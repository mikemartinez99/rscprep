# flagQC()'s actual flagging: where the thresholds land, which cells they
# catch, and what the summary table records.
#
# Fixtures pin the statistics rather than approximate them. metric_values()
# has a median of 6 and a MAD of exactly 2 * 1.4826, and appending one extreme
# value shifts neither, so the expected thresholds are arithmetic:
# 6 +/- 3 * 2 * 1.4826 at the default nmads.
#
# Regexes use \\s+ where a message could wrap, since cli hard-wraps at the
# console width.

# ---- thresholds ----------------------------------------------------------

test_that("thresholds are median +/- nmads * MAD on the raw scale", {
    obj <- make_metric_seurat(m = metric_values(100))

    qc <- flagQC(obj, c(m = "both"), log = FALSE)@misc$flagQC

    expect_equal(qc$center, 6)
    expect_equal(qc$lower, 6 - 3 * 2 * MAD_CONST)
    expect_equal(qc$upper, 6 + 3 * 2 * MAD_CONST)
})

test_that("nmads scales the thresholds", {
    obj <- make_metric_seurat(m = metric_values(100))

    qc <- flagQC(obj, c(m = "both"), nmads = 5, log = FALSE)@misc$flagQC

    expect_equal(qc$upper, 6 + 5 * 2 * MAD_CONST)
})

test_that("the flagged cells are exactly those outside the reported thresholds", {
    values <- metric_values(100)
    obj    <- make_metric_seurat(m = values)

    res <- flagQC(obj, c(m = "both"), log = FALSE)
    qc  <- res@misc$flagQC

    expect_equal(unname(res$m_3_both_outliers),
                 values < qc$lower | values > qc$upper)
})

test_that("counts are split between the two tails", {
    obj <- make_metric_seurat(m = metric_values(100))

    qc <- flagQC(obj, c(m = "both"), log = FALSE)@misc$flagQC

    expect_equal(qc$n_lower, 0L)
    expect_equal(qc$n_upper, 1L)
    expect_equal(qc$n_flagged, 1L)
})

test_that("the two tail counts sum to the total flagged", {
    obj <- make_metric_seurat(m = c(metric_values(100), -100))

    qc <- flagQC(obj, c(m = "both"), log = FALSE)@misc$flagQC

    expect_equal(qc$n_lower + qc$n_upper, qc$n_flagged)
    expect_equal(qc$n_flagged, 2L)
})

test_that("a distribution with no outliers flags nothing", {
    obj <- make_metric_seurat(m = metric_values())

    res <- flagQC(obj, c(m = "both"), log = FALSE)

    expect_equal(sum(res$m_3_both_outliers), 0L)
    expect_false(any(is.na(res$m_3_both_outliers)))
})

# ---- bounds --------------------------------------------------------------

test_that("bound = 'lower' flags only the low tail", {
    obj <- make_metric_seurat(m = c(metric_values(100), -100))

    res <- flagQC(obj, c(m = "lower"), log = FALSE)
    qc  <- res@misc$flagQC

    expect_equal(qc$n_flagged, 1L)
    expect_true(res$m_3_lower_outliers[22])   # the -100
    expect_false(res$m_3_lower_outliers[21])  # the 100
})

test_that("bound = 'higher' flags only the high tail", {
    obj <- make_metric_seurat(m = c(metric_values(100), -100))

    res <- flagQC(obj, c(m = "higher"), log = FALSE)

    expect_equal(res@misc$flagQC$n_flagged, 1L)
    expect_true(res$m_3_higher_outliers[21])
    expect_false(res$m_3_higher_outliers[22])
})

test_that("an unenforced tail is NA, not an unapplied number", {
    obj <- make_metric_seurat(m = metric_values(100))

    lo <- flagQC(obj, c(m = "lower"), log = FALSE)@misc$flagQC
    hi <- flagQC(obj, c(m = "higher"), log = FALSE)@misc$flagQC

    expect_true(is.na(lo$upper))
    expect_true(is.na(lo$n_upper))
    expect_false(is.na(lo$lower))

    expect_true(is.na(hi$lower))
    expect_true(is.na(hi$n_lower))
    expect_false(is.na(hi$upper))
})

test_that("the bound is recorded in the column name and the table", {
    obj <- make_metric_seurat(m = metric_values(100))

    res <- flagQC(obj, c(m = "higher"), nmads = 2, log = FALSE)

    expect_equal(res@misc$flagQC$column, "m_2_higher_outliers")
    expect_equal(res@misc$flagQC$bound, "higher")
    expect_equal(res@misc$flagQC$nmads, 2)
})

# ---- the log transform ---------------------------------------------------

test_that("thresholds are estimated on the log scale and reported on the raw one", {
    obj <- make_metric_seurat(m = log_metric_values(1024))

    qc <- flagQC(obj, c(m = "both"), log = TRUE)@misc$flagQC

    expect_equal(qc$scale, "log2")
    expect_equal(qc$center, 2^3)
    expect_equal(qc$upper, 2^(3 + 3 * MAD_CONST))
    expect_equal(qc$lower, 2^(3 - 3 * MAD_CONST))
})

test_that("log flags exactly the cells outside the back-transformed thresholds", {
    values <- log_metric_values(1024)
    obj    <- make_metric_seurat(m = values)

    res <- flagQC(obj, c(m = "both"), log = TRUE)
    qc  <- res@misc$flagQC

    expect_equal(unname(res$m_3_both_outliers),
                 values < qc$lower | values > qc$upper)
    expect_equal(qc$n_upper, 1L)
})

test_that("log = FALSE is recorded as the identity scale", {
    obj <- make_metric_seurat(m = metric_values(100))

    expect_equal(flagQC(obj, "m", log = FALSE)@misc$flagQC$scale, "identity")
})

test_that("a metric containing zeroes is shifted and shifted back", {
    # log2(x + 1) on the way in must be undone on the way out, or every
    # reported threshold is off by one
    values <- c(0, 0, metric_values(100))
    obj    <- make_metric_seurat(m = values)

    res <- flagQC(obj, c(m = "both"), log = TRUE)
    qc  <- res@misc$flagQC

    expect_equal(qc$center, stats::median(values))
    expect_equal(unname(res$m_3_both_outliers),
                 values < qc$lower | values > qc$upper)
})

test_that("the log and identity scales give different thresholds", {
    obj <- make_metric_seurat(m = log_metric_values(1024))

    logged <- flagQC(obj, "m", log = TRUE)@misc$flagQC
    raw    <- flagQC(obj, "m", log = FALSE)@misc$flagQC

    expect_false(isTRUE(all.equal(logged$upper, raw$upper)))
})

# ---- missing values ------------------------------------------------------

test_that("cells with a missing metric are left NA, not FALSE", {
    values <- metric_values(100)
    values[c(1, 2)] <- NA
    obj <- make_metric_seurat(m = values)

    res <- flagQC(obj, c(m = "both"), log = FALSE)

    expect_true(all(is.na(res$m_3_both_outliers[1:2])))
    expect_false(any(is.na(res$m_3_both_outliers[-(1:2)])))
})

test_that("missing values are counted in the table", {
    values <- metric_values(100)
    values[c(1, 2)] <- NA

    qc <- flagQC(make_metric_seurat(m = values), c(m = "both"),
                 log = FALSE)@misc$flagQC

    expect_equal(qc$n_missing, 2L)
    expect_equal(qc$n_cells, 21L)
})

test_that("missing values do not shift the thresholds", {
    values     <- metric_values(100)
    with_na    <- values
    with_na[1] <- NA

    a <- flagQC(make_metric_seurat(m = values), "m", log = FALSE)@misc$flagQC
    b <- flagQC(make_metric_seurat(m = with_na), "m", log = FALSE)@misc$flagQC

    expect_equal(a$center, b$center)
})

# ---- degenerate metrics --------------------------------------------------

test_that("a zero-MAD metric warns and writes no column", {
    # over half the cells share one value, so both thresholds land on it
    obj <- make_metric_seurat(m = rep(5, 20))

    expect_warning(flagQC(obj, "m", log = FALSE), "no threshold")
    res <- suppressWarnings(flagQC(obj, "m", log = FALSE))

    expect_false("m_3_both_outliers" %in% colnames(res@meta.data))
})

test_that("a skipped metric still appears in the table, with NA thresholds", {
    obj <- make_metric_seurat(m = rep(5, 20))

    qc <- suppressWarnings(flagQC(obj, "m", log = FALSE))@misc$flagQC

    expect_equal(nrow(qc), 1)
    expect_equal(qc$metric, "m")
    expect_true(is.na(qc$lower))
    expect_true(is.na(qc$upper))
    expect_true(is.na(qc$column))
})

test_that("one degenerate metric does not stop the others", {
    obj <- make_metric_seurat(flat = rep(5, 21), m = metric_values(100))

    res <- suppressWarnings(flagQC(obj, c("flat", "m"), log = FALSE))

    expect_false("flat_3_both_outliers" %in% colnames(res@meta.data))
    expect_true("m_3_both_outliers" %in% colnames(res@meta.data))
    expect_equal(nrow(res@misc$flagQC), 2)
})

# ---- the summary table ---------------------------------------------------

test_that("the summary table has one row per metric", {
    obj <- make_metric_seurat(a = metric_values(100), b = metric_values(100))

    qc <- flagQC(obj, c("a", "b"), log = FALSE)@misc$flagQC

    expect_s3_class(qc, "data.frame")
    expect_equal(nrow(qc), 2)
    expect_equal(qc$metric, c("a", "b"))
})

test_that("the summary table has the documented columns", {
    obj <- make_metric_seurat(m = metric_values(100))

    qc <- flagQC(obj, "m", log = FALSE)@misc$flagQC

    expect_named(qc, c("metric", "bound", "nmads", "scale", "group_by", "group",
                       "center", "lower", "upper", "n_lower", "n_upper",
                       "n_flagged", "pct_flagged", "n_missing", "n_cells",
                       "column"))
})

test_that("group columns are NA when thresholds are global", {
    obj <- make_metric_seurat(m = metric_values(100))

    qc <- flagQC(obj, "m", log = FALSE)@misc$flagQC

    expect_true(is.na(qc$group_by))
    expect_true(is.na(qc$group))
})

test_that("pct_flagged is the flagged share of the cells", {
    obj <- make_metric_seurat(m = metric_values(100))

    qc <- flagQC(obj, "m", log = FALSE)@misc$flagQC

    expect_equal(qc$pct_flagged, 100 * qc$n_flagged / qc$n_cells)
})

test_that("repeated calls accumulate into one table", {
    obj <- make_metric_seurat(a = metric_values(100), b = metric_values(100))

    res <- flagQC(obj, "a", log = FALSE)
    res <- flagQC(res, "b", log = FALSE)

    expect_equal(nrow(res@misc$flagQC), 2)
    expect_equal(res@misc$flagQC$metric, c("a", "b"))
})

test_that("re-flagging the same column replaces its rows rather than duplicating", {
    obj <- make_metric_seurat(m = metric_values(100))

    res <- flagQC(obj, "m", log = FALSE)
    res <- flagQC(res, "m", log = FALSE)

    expect_equal(nrow(res@misc$flagQC), 1)
})

test_that("every table row corresponds to a column on the object", {
    obj <- make_metric_seurat(a = metric_values(100), b = metric_values(100))

    res <- flagQC(obj, c(a = "lower", b = "higher"), log = FALSE)
    qc  <- res@misc$flagQC

    expect_true(all(stats::na.omit(qc$column) %in% colnames(res@meta.data)))
})

# ---- what is returned ----------------------------------------------------

test_that("a Seurat object is returned with nothing filtered", {
    obj <- make_metric_seurat(m = metric_values(100))

    res <- flagQC(obj, "m", log = FALSE)

    expect_s4_class(res, "Seurat")
    expect_equal(ncol(res), ncol(obj))
    expect_equal(colnames(res), colnames(obj))
})

test_that("the flag column is logical", {
    obj <- make_metric_seurat(m = metric_values(100))

    res <- flagQC(obj, "m", log = FALSE)

    expect_type(res$m_3_both_outliers, "logical")
})

test_that("existing metadata is untouched", {
    obj <- make_metric_seurat(m = metric_values(100))

    res <- flagQC(obj, "m", log = FALSE)

    expect_equal(res$m, obj$m)
    expect_equal(res$Ident, obj$Ident)
})
