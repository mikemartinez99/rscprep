# flagQC(group_by=): thresholds estimated within groups rather than across the
# whole object.
#
# The case that motivates the argument is two samples at different sequencing
# depths. A global threshold sits between them and describes neither; a
# per-group threshold describes each. Several tests below are built on exactly
# that shape, using metric_values() scaled per group so the expected cutoffs
# stay hand-computable within each one.
#
# Regexes use \\s+ where a message could wrap, since cli hard-wraps at the
# console width.

#' Two groups an order of magnitude apart, each with its own high outlier.
#'
#' Group A is metric_values() as-is (median 6); group B is the same values
#' times 100 (median 600). Each carries one extreme, so per-group flagging
#' catches two cells and a global threshold catches neither tail cleanly.
two_scale_seurat <- function() {
    make_metric_seurat(
        m     = c(metric_values(100), metric_values(100) * 100),
        ident = rep(c("A", "B"), each = 21)
    )
}

# ---- one row per group ---------------------------------------------------

test_that("grouping produces one table row per group", {
    res <- flagQC(two_scale_seurat(), c(m = "both"),
                  group_by = "Ident", log = FALSE)
    qc  <- res@misc$flagQC

    expect_equal(nrow(qc), 2)
    expect_equal(qc$group, c("A", "B"))
    expect_true(all(qc$group_by == "Ident"))
})

test_that("each group gets its own thresholds", {
    qc <- flagQC(two_scale_seurat(), c(m = "both"),
                 group_by = "Ident", log = FALSE)@misc$flagQC

    expect_equal(qc$center, c(6, 600))
    expect_equal(qc$upper, c(6 + 3 * 2 * MAD_CONST, 600 + 3 * 200 * MAD_CONST))
})

test_that("group cell counts sum to the object's", {
    res <- flagQC(two_scale_seurat(), "m", group_by = "Ident", log = FALSE)

    expect_equal(sum(res@misc$flagQC$n_cells), ncol(res))
})

test_that("pct_flagged is within-group, not object-wide", {
    qc <- flagQC(two_scale_seurat(), c(m = "higher"),
                 group_by = "Ident", log = FALSE)@misc$flagQC

    expect_true(all(qc$n_cells == 21))
    expect_equal(qc$pct_flagged, 100 * qc$n_flagged / 21)
})

# ---- the case grouping exists for ----------------------------------------

test_that("a global threshold describes neither group when the scales differ", {
    obj <- two_scale_seurat()

    global  <- flagQC(obj, c(m = "higher"), log = FALSE)
    grouped <- flagQC(obj, c(m = "higher"), group_by = "Ident", log = FALSE)

    # grouped: each group's own extreme, and nothing else
    expect_equal(unname(which(grouped$m_3_higher_outliers)), c(21L, 42L))

    # global: the single threshold lands inside group B's ordinary range, so
    # it misses group A's outlier entirely and flags four healthy B cells
    # instead -- wrong in both directions at once
    expect_false(global$m_3_higher_outliers[[21]])
    expect_equal(sum(global$m_3_higher_outliers[22:41]), 4L)
})

test_that("a cell is judged against its own group", {
    obj <- two_scale_seurat()

    res <- flagQC(obj, c(m = "higher"), group_by = "Ident", log = FALSE)
    qc  <- res@misc$flagQC

    flagged <- which(res$m_3_higher_outliers)
    expect_equal(unname(flagged), c(21L, 42L))   # each group's extreme
    expect_equal(qc$n_upper, c(1L, 1L))
})

test_that("thresholds match what the same values give ungrouped", {
    # a group's threshold should be exactly what flagQC would produce if that
    # group were the whole object
    obj <- two_scale_seurat()

    grouped <- flagQC(obj, c(m = "both"), group_by = "Ident",
                      log = FALSE)@misc$flagQC
    alone   <- flagQC(make_metric_seurat(m = metric_values(100)), c(m = "both"),
                      log = FALSE)@misc$flagQC

    expect_equal(grouped$lower[1], alone$lower)
    expect_equal(grouped$upper[1], alone$upper)
})

# ---- what counts as a group ----------------------------------------------

test_that("a single-level grouping column works", {
    obj <- make_metric_seurat(m = metric_values(100), ident = "only")

    qc <- flagQC(obj, "m", group_by = "Ident", log = FALSE)@misc$flagQC

    expect_equal(nrow(qc), 1)
    expect_equal(qc$group, "only")
})

test_that("a single-level grouping gives the same thresholds as no grouping", {
    obj <- make_metric_seurat(m = metric_values(100), ident = "only")

    grouped <- flagQC(obj, "m", group_by = "Ident", log = FALSE)@misc$flagQC
    global  <- flagQC(obj, "m", log = FALSE)@misc$flagQC

    expect_equal(grouped$lower, global$lower)
    expect_equal(grouped$upper, global$upper)
})

test_that("integer group codes are accepted", {
    obj <- two_scale_seurat()
    obj$code <- rep(c(1L, 2L), each = 21)

    qc <- flagQC(obj, "m", group_by = "code", log = FALSE)@misc$flagQC

    expect_equal(nrow(qc), 2)
    expect_equal(qc$group, c("1", "2"))
})

test_that("factor levels with no cells are dropped", {
    obj <- two_scale_seurat()
    obj$fct <- factor(rep(c("A", "B"), each = 21), levels = c("A", "B", "unused"))

    qc <- flagQC(obj, "m", group_by = "fct", log = FALSE)@misc$flagQC

    expect_equal(nrow(qc), 2)
    expect_false("unused" %in% qc$group)
})

test_that("any metadata column can be grouped on, not just Ident", {
    obj <- two_scale_seurat()
    obj$batch <- rep(c("b1", "b2"), times = 21)

    qc <- flagQC(obj, "m", group_by = "batch", log = FALSE)@misc$flagQC

    expect_equal(qc$group, c("b1", "b2"))
    expect_true(all(qc$group_by == "batch"))
})

# ---- cells with no group -------------------------------------------------

test_that("cells with a missing group value are warned about", {
    obj <- two_scale_seurat()
    obj$grp <- c(rep("A", 20), NA, rep("B", 21))

    expect_warning(flagQC(obj, "m", group_by = "grp", log = FALSE),
                   "no grp value")
})

test_that("cells with a missing group value are left NA, not FALSE", {
    obj <- two_scale_seurat()
    obj$grp <- c(rep("A", 20), NA, rep("B", 21))

    res <- suppressWarnings(flagQC(obj, "m", group_by = "grp", log = FALSE))

    expect_true(is.na(res$m_3_both_outliers[21]))
    expect_false(any(is.na(res$m_3_both_outliers[-21])))
})

test_that("ungrouped cells are excluded from every group's estimate", {
    obj <- two_scale_seurat()
    obj$grp <- c(rep("A", 21), rep(NA, 21))

    qc <- suppressWarnings(
        flagQC(obj, "m", group_by = "grp", log = FALSE)
    )@misc$flagQC

    expect_equal(nrow(qc), 1)
    expect_equal(qc$n_cells, 21L)
})

# ---- groups too small or too flat to estimate from -----------------------

test_that("a group with fewer than three cells is skipped", {
    obj <- make_metric_seurat(
        m     = c(metric_values(100), 1, 2),
        ident = c(rep("big", 21), "tiny", "tiny")
    )

    expect_warning(flagQC(obj, "m", group_by = "Ident", log = FALSE),
                   "fewer than 3 cells")
})

test_that("a skipped group's cells are left NA, not FALSE", {
    obj <- make_metric_seurat(
        m     = c(metric_values(100), 1, 2),
        ident = c(rep("big", 21), "tiny", "tiny")
    )

    res <- suppressWarnings(flagQC(obj, "m", group_by = "Ident", log = FALSE))

    expect_true(all(is.na(res$m_3_both_outliers[22:23])))
    expect_false(any(is.na(res$m_3_both_outliers[1:21])))
})

test_that("a skipped group still appears in the table", {
    obj <- make_metric_seurat(
        m     = c(metric_values(100), 1, 2),
        ident = c(rep("big", 21), "tiny", "tiny")
    )

    qc <- suppressWarnings(
        flagQC(obj, "m", group_by = "Ident", log = FALSE)
    )@misc$flagQC

    expect_equal(nrow(qc), 2)
    tiny <- qc[qc$group == "tiny", ]
    expect_true(is.na(tiny$lower))
    expect_true(is.na(tiny$column))
    expect_equal(tiny$n_cells, 2L)
})

test_that("a zero-MAD group is skipped while its neighbours are not", {
    obj <- make_metric_seurat(
        m     = c(metric_values(100), rep(5, 21)),
        ident = rep(c("varied", "flat"), each = 21)
    )

    expect_warning(flagQC(obj, "m", group_by = "Ident", log = FALSE),
                   "zero MAD")

    res <- suppressWarnings(flagQC(obj, "m", group_by = "Ident", log = FALSE))
    expect_true(all(is.na(res$m_3_both_outliers[22:42])))
    expect_equal(sum(res$m_3_both_outliers[1:21]), 1L)
})

test_that("skipped groups are named in one warning, not one warning each", {
    obj <- make_metric_seurat(
        m     = c(metric_values(100), 1, 2, 3, 4),
        ident = c(rep("big", 21), "t1", "t1", "t2", "t2")
    )

    w <- testthat::capture_warnings(flagQC(obj, "m", group_by = "Ident",
                                           log = FALSE))

    expect_length(w, 1)
    expect_match(w, "t1")
    expect_match(w, "t2")
    expect_match(w, "2 of 3")
})

test_that("no column is written when every group is skipped", {
    obj <- make_metric_seurat(m = c(1, 2, 3, 4),
                              ident = c("a", "a", "b", "b"))

    res <- suppressWarnings(flagQC(obj, "m", group_by = "Ident", log = FALSE))

    expect_false("m_3_both_outliers" %in% colnames(res@meta.data))
    expect_true(all(is.na(res@misc$flagQC$column)))
})

# ---- interaction with the rest of the arguments --------------------------

test_that("grouping composes with per-metric bounds", {
    obj <- two_scale_seurat()
    obj$other <- c(metric_values(100), metric_values(100) * 100)

    qc <- flagQC(obj, c(m = "lower", other = "higher"),
                 group_by = "Ident", log = FALSE)@misc$flagQC

    expect_equal(nrow(qc), 4)
    expect_true(all(is.na(qc$upper[qc$metric == "m"])))
    expect_true(all(is.na(qc$lower[qc$metric == "other"])))
})

test_that("grouping composes with the log transform", {
    obj <- make_metric_seurat(
        m     = c(log_metric_values(1024), log_metric_values(1024) * 8),
        ident = rep(c("A", "B"), each = 21)
    )

    qc <- flagQC(obj, c(m = "both"), group_by = "Ident", log = TRUE)@misc$flagQC

    expect_equal(qc$scale, c("log2", "log2"))
    expect_equal(qc$center, c(2^3, 2^6))
    expect_equal(qc$upper, c(2^(3 + 3 * MAD_CONST), 2^(6 + 3 * MAD_CONST)))
})

test_that("re-flagging with different grouping replaces the table rows", {
    obj <- two_scale_seurat()

    res <- flagQC(obj, "m", log = FALSE)
    expect_true(is.na(res@misc$flagQC$group_by))

    res <- flagQC(res, "m", group_by = "Ident", log = FALSE)

    expect_equal(nrow(res@misc$flagQC), 2)
    expect_true(all(res@misc$flagQC$group_by == "Ident"))
})
