# pointAesthetics() and the figure it feeds.
#
# A barcode rank curve can be a few hundred points or a few thousand depending
# on how many distinct totals the raw matrix holds, so point size and opacity
# are derived from that count rather than fixed. What is asserted here is the
# shape of that relationship -- denser curves get smaller, more translucent
# points -- and that the values actually reach the rendered layer.
#
# Same fixture strategy as the other reader tests: real directory layouts with
# placeholder H5 files and Seurat::Read10X_h5() mocked.

# ---- pointAesthetics -----------------------------------------------------

test_that("both aesthetics are returned", {
    p <- pointAesthetics(1000)

    expect_named(p, c("size", "alpha"))
    expect_type(p$size, "double")
    expect_type(p$alpha, "double")
})

test_that("denser curves get smaller, more translucent points", {
    sparse <- pointAesthetics(200)
    dense  <- pointAesthetics(5000)

    expect_gt(sparse$size, dense$size)
    expect_gt(sparse$alpha, dense$alpha)
})

test_that("size and alpha are monotone in the point count", {
    n     <- c(1, 10, 100, 1000, 10000, 1e6)
    size  <- vapply(n, function(i) pointAesthetics(i)$size, numeric(1))
    alpha <- vapply(n, function(i) pointAesthetics(i)$alpha, numeric(1))

    expect_false(is.unsorted(rev(size)))
    expect_false(is.unsorted(rev(alpha)))
})

test_that("both are clamped outside the anchor range", {
    # below 100 points and above 10,000 the values plateau rather than running
    # off to invisible or opaque
    expect_equal(pointAesthetics(1), pointAesthetics(100))
    expect_equal(pointAesthetics(10000), pointAesthetics(1e9))
})

test_that("alpha stays a usable opacity at every size", {
    for (n in c(0, 1, 100, 1000, 10000, 1e9)) {
        a <- pointAesthetics(n)$alpha
        expect_true(a > 0 && a <= 1)
    }
    expect_true(all(vapply(c(0, 1, 1e9),
                           function(n) pointAesthetics(n)$size > 0, TRUE)))
})

test_that("degenerate point counts do not error", {
    # nrow() should never hand over these, but an empty curve must not take
    # the whole figure down with it
    expect_type(pointAesthetics(0)$size, "double")
    expect_type(pointAesthetics(NA)$size, "double")
    expect_type(pointAesthetics(numeric(0))$size, "double")
})

# ---- the figure it feeds -------------------------------------------------

test_that("the plot builds, which is what a missing pointAesthetics broke", {
    root <- local_root()
    local_mock_read10x(counts = make_raw_counts())
    outs <- make_count_run(root, "Sample_1")
    add_raw_h5(outs, "flat")

    res <- plotBarcodeRanks(path = outs, multi = FALSE)

    expect_named(res, "Sample_1")
    expect_s3_class(res$Sample_1$plot, "ggplot")
    # constructing a ggplot defers evaluation; building it is what actually
    # resolves the layers
    expect_no_error(ggplot2::ggplot_build(res$Sample_1$plot))
})

test_that("the computed aesthetics reach the rendered point layer", {
    root <- local_root()
    local_mock_read10x(counts = make_raw_counts())
    outs <- make_count_run(root, "Sample_1")
    add_raw_h5(outs, "flat")

    res   <- plotBarcodeRanks(path = outs, multi = FALSE)
    built <- ggplot2::ggplot_build(res$Sample_1$plot)
    point <- pointAesthetics(nrow(res$Sample_1$data))

    expect_equal(unique(built$data[[1]]$size), point$size)
    expect_equal(unique(built$data[[1]]$alpha), point$alpha)
})

# ---- DropletUtils versions -----------------------------------------------
#
# barcodeRanks() stopped returning a `fitted` column in DropletUtils 1.30.
# The real function only ever shows one shape -- whichever version is
# installed -- so both are built by hand here.

fake_ranked <- function(fitted = TRUE) {
    skip_if_not_installed("S4Vectors")
    totals <- c(5000, 3000, 800, 40, 10, 0)
    cols   <- list(rank = seq_along(totals), total = totals)
    if (fitted) cols$fitted <- c(NA, 3100, 790, 42, NA, NA)
    out <- do.call(S4Vectors::DataFrame, cols)
    rownames(out) <- paste0("BC", seq_along(totals))
    S4Vectors::metadata(out) <- list(knee = 3000, inflection = 800)
    out
}

test_that("a barcodeRanks() result with no fitted column is summarised (DropletUtils >= 1.30)", {
    res <- suppressMessages(summariseBarcodeRanks(fake_ranked(fitted = FALSE), "S1", "x.h5"))

    expect_equal(nrow(res$data), 5)          # the zero-total barcode is dropped
    expect_true(all(is.na(res$data$fitted)))
    expect_equal(res$n_above_inflection, 3)
})

test_that("a barcodeRanks() result with a fitted column keeps it (DropletUtils < 1.30)", {
    res <- suppressMessages(summariseBarcodeRanks(fake_ranked(fitted = TRUE), "S1", "x.h5"))

    expect_equal(res$data$fitted, c(NA, 3100, 790, 42, NA))
})

test_that("the fitted line is drawn only when there is a spline to draw", {
    with_fit    <- suppressMessages(summariseBarcodeRanks(fake_ranked(TRUE), "S1", "x.h5"))
    without_fit <- suppressMessages(summariseBarcodeRanks(fake_ranked(FALSE), "S1", "x.h5"))

    n_lines <- function(p) sum(vapply(p$layers, function(l) inherits(l$geom, "GeomLine"), logical(1)))
    expect_equal(n_lines(barcodeRankPlot(with_fit, "S1")), 1)
    expect_equal(n_lines(barcodeRankPlot(without_fit, "S1")), 0)
})
