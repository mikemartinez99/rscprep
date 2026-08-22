# Read10X_h5() is mocked to hand back a synthetic raw matrix with a real knee
# in it, but DropletUtils::barcodeRanks() is NOT mocked -- it runs for real, so
# these tests also verify the matrix is handed to it in a shape it accepts.

skip_if_no_droplet_utils <- function() {
    testthat::skip_if_not_installed("DropletUtils")
    testthat::skip_if_not_installed("S4Vectors")
}

#' A `cellranger count` run whose raw H5 reads back as `counts`.
local_rank_run <- function(root, sample = "Sample_1",
                           counts = make_raw_counts(), .env = parent.frame()) {
    local_mock_read10x(counts = counts, .env = .env)
    outs <- make_count_run(root, sample)
    add_raw_h5(outs, "flat")
    outs
}

test_that("one run returns a curve, a plot and both cutoffs", {
    skip_if_no_droplet_utils()
    root <- local_root()
    outs <- local_rank_run(root)

    res <- plotBarcodeRanks(path = outs, multi = FALSE)

    expect_named(res, "Sample_1")
    expect_named(res$Sample_1,
                 c("data", "knee", "inflection", "n_barcodes",
                   "n_above_inflection", "h5", "plot"),
                 ignore.order = TRUE)
    expect_s3_class(res$Sample_1$plot, "ggplot")
    expect_true(res$Sample_1$knee >= res$Sample_1$inflection)
    expect_equal(res$Sample_1$h5, file.path(outs, "raw_feature_bc_matrix.h5"))
})

test_that("the knee separates the cell population from the ambient tail", {
    skip_if_no_droplet_utils()
    root <- local_root()
    # 300 real cells at ~40 UMIs/gene over 40 genes, then 3000 ambient barcodes
    outs <- local_rank_run(root, counts = make_raw_counts(n_cells = 300,
                                                          n_empty = 3000))

    res <- plotBarcodeRanks(path = outs, multi = FALSE)$Sample_1

    expect_equal(res$n_barcodes, 3300)
    # the cutoff should land near the 300 planted cells, not out in the tail
    expect_gt(res$n_above_inflection, 200)
    expect_lt(res$n_above_inflection, 500)
})

test_that("the returned curve is the plotted curve, not the full barcode table", {
    skip_if_no_droplet_utils()
    root <- local_root()
    outs <- local_rank_run(root)

    curve <- plotBarcodeRanks(path = outs, multi = FALSE)$Sample_1$data

    expect_named(curve, c("barcode", "rank", "total", "fitted"))
    # one row per distinct total, so far fewer rows than barcodes
    expect_lt(nrow(curve), 3300)
    expect_equal(anyDuplicated(curve$rank), 0)
    # nothing at zero -- a log axis cannot show it
    expect_true(all(curve$total > 0))
    expect_false(is.unsorted(curve$rank))
})

test_that("`lower` is passed through to barcodeRanks()", {
    skip_if_no_droplet_utils()
    root <- local_root()
    outs <- local_rank_run(root)

    default <- plotBarcodeRanks(path = outs, multi = FALSE)$Sample_1
    # high enough to bite: this fixture's ambient tail all sits below 100, so a
    # modest bump would exclude nothing and change nothing
    raised  <- plotBarcodeRanks(path = outs, multi = FALSE,
                                lower = 1000)$Sample_1

    # a different ambient assumption moves the fit, hence the cutoffs
    expect_false(isTRUE(all.equal(default$inflection, raised$inflection)))
    expect_lt(sum(!is.na(raised$data$fitted)), sum(!is.na(default$data$fitted)))
    # ...but not the barcodes themselves: `lower` is not a filter
    expect_equal(default$n_barcodes, raised$n_barcodes)
})

test_that("extra arguments reach barcodeRanks()", {
    skip_if_no_droplet_utils()
    root <- local_root()
    outs <- local_rank_run(root)

    default <- plotBarcodeRanks(path = outs, multi = FALSE)$Sample_1
    stiffer <- plotBarcodeRanks(path = outs, multi = FALSE, df = 5)$Sample_1

    expect_false(isTRUE(all.equal(default$data$fitted, stiffer$data$fitted)))
})

test_that("outdir writes one CSV and one PNG per run", {
    skip_if_no_droplet_utils()
    root <- local_root()
    outs <- local_rank_run(root)
    dest <- file.path(root, "qc", "ranks")   # does not exist yet

    res <- plotBarcodeRanks(path = outs, multi = FALSE, outdir = dest,
                            width = 5, height = 4, dpi = 72)

    expect_setequal(list.files(dest),
                    c("Sample_1_barcode_ranks.csv", "Sample_1_barcode_ranks.png"))

    written <- utils::read.csv(file.path(dest, "Sample_1_barcode_ranks.csv"))
    expect_equal(nrow(written), nrow(res$Sample_1$data))
    expect_equal(written$total, res$Sample_1$data$total)
})

test_that("nothing is written when outdir is NULL", {
    skip_if_no_droplet_utils()
    root <- local_root()
    outs <- local_rank_run(root)
    before <- list.files(root, recursive = TRUE)

    plotBarcodeRanks(path = outs, multi = FALSE)

    expect_equal(list.files(root, recursive = TRUE), before)
})

test_that("run names are sanitised for use as filenames", {
    skip_if_no_droplet_utils()
    root <- local_root()
    outs <- local_rank_run(root, sample = "Pool 1/rerun")
    dest <- file.path(root, "out")

    plotBarcodeRanks(path = outs, multi = FALSE, outdir = dest)

    expect_setequal(list.files(dest),
                    c("rerun_barcode_ranks.csv", "rerun_barcode_ranks.png"))
})

test_that("a missing raw matrix is fatal in path mode and skipped in manifest mode", {
    skip_if_no_droplet_utils()
    root    <- local_root()
    present <- local_rank_run(root, "Sample_1")
    absent  <- make_count_run(root, "Sample_2")   # filtered only, no raw

    expect_error(plotBarcodeRanks(path = absent, multi = FALSE),
                 "Raw matrix for")

    manifest <- data.frame(Path  = c(absent, present),
                           multi = c(FALSE, FALSE))
    expect_warning(res <- plotBarcodeRanks(manifest = manifest), "Raw matrix for")
    expect_named(res, "Sample_1")
})

test_that("a run barcodeRanks() cannot fit is skipped, not fatal, in manifest mode", {
    skip_if_no_droplet_utils()
    root <- local_root()

    # too few distinct totals for a spline fit
    degenerate <- make_counts(n_genes = 5, n_cells = 4, extras = FALSE)
    local_mock_read10x(counts = degenerate)

    bad <- make_count_run(root, "Sample_Bad")
    add_raw_h5(bad, "flat")

    expect_warning(
        res <- plotBarcodeRanks(manifest = data.frame(Path = bad,
                                                      multi = FALSE)),
        "barcode ranking failed"
    )
    expect_length(res, 0)
})

test_that("the plot is a plain ggplot with log-log axes and both cutoff lines", {
    skip_if_no_droplet_utils()
    root <- local_root()
    outs <- local_rank_run(root)

    res <- plotBarcodeRanks(path = outs, multi = FALSE)$Sample_1
    p   <- res$plot

    # points, fitted line, two dashed rules, two labels
    expect_length(p$layers, 6)
    expect_equal(p$labels$title, "Sample_1")
    expect_equal(p$labels$x, "Barcode rank")

    # both annotations survive the log transform -- x = -Inf would not
    built <- ggplot2::ggplot_build(p)
    labels <- unlist(lapply(built$data, function(d) d$label))
    expect_true(any(grepl("^Knee: ", labels)))
    expect_true(any(grepl("^Inflection: ", labels)))
    expect_true(all(is.finite(unlist(lapply(built$data, function(d) d$x)))))
})

test_that("argument validation matches seuratFrom10x()", {
    skip_if_no_droplet_utils()

    expect_error(plotBarcodeRanks(), "Provide either `manifest` or `path`")
    expect_error(plotBarcodeRanks(path = "/some/outs"),
                 "you must also specify `multi`")
    expect_error(plotBarcodeRanks(manifest = data.frame(other = 1)),
                 "missing required columns: Path and multi")
    expect_error(plotBarcodeRanks(path = "/x", multi = FALSE, lower = -1),
                 "'lower' must be a single non-negative number")
    expect_error(plotBarcodeRanks(path = "/x", multi = FALSE, outdir = ""),
                 "'outdir' must be a single directory path")
})
