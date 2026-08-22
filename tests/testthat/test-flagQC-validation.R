# flagQC()'s argument surface. Everything is checked before any metric is
# touched, so a malformed call fails without half-writing columns.
#
# Regexes use \\s+ where a message could wrap, since cli hard-wraps at the
# console width.

# ---- the object ----------------------------------------------------------

test_that("a list of objects is rejected with a usable hint", {
    obj_list <- make_seurat_list()

    expect_error(flagQC(obj_list), "is a list, not a single Seurat object")
    expect_error(flagQC(obj_list), "lapply\\(obj_list, flagQC")
})

test_that("non-Seurat input is rejected", {
    expect_error(flagQC(42), "must be a Seurat object")
})

test_that("errors name the argument as the caller wrote it", {
    my_object <- 42
    expect_error(flagQC(my_object), "my_object")
})

# ---- nmads ---------------------------------------------------------------

test_that("nmads must be a single positive number", {
    obj <- make_metric_seurat(m = metric_values())

    expect_error(flagQC(obj, "m", nmads = 0), "single positive number")
    expect_error(flagQC(obj, "m", nmads = -1), "single positive number")
    expect_error(flagQC(obj, "m", nmads = NA), "single positive number")
    expect_error(flagQC(obj, "m", nmads = "3"), "single positive number")
    expect_error(flagQC(obj, "m", nmads = c(2, 3)), "single positive number")
})

test_that("a fractional nmads is accepted", {
    obj <- make_metric_seurat(m = metric_values(100))
    expect_no_error(flagQC(obj, "m", nmads = 2.5, log = FALSE))
})

# ---- log -----------------------------------------------------------------

test_that("log must be a single non-NA logical", {
    obj <- make_metric_seurat(m = metric_values())

    expect_error(flagQC(obj, "m", log = "yes"), "TRUE.*FALSE")
    expect_error(flagQC(obj, "m", log = NA), "TRUE.*FALSE")
    expect_error(flagQC(obj, "m", log = c(TRUE, FALSE)), "TRUE.*FALSE")
})

# ---- metrics -------------------------------------------------------------

test_that("an empty metrics argument is an error", {
    obj <- make_metric_seurat(m = metric_values())
    expect_error(flagQC(obj, character(0)), "is empty")
})

test_that("a non-character metrics argument is an error", {
    obj <- make_metric_seurat(m = metric_values())
    expect_error(flagQC(obj, 1:2), "must be a character vector")
})

test_that("absent metadata columns are named", {
    obj <- make_metric_seurat(m = metric_values())

    expect_error(flagQC(obj, c("m", "nope")), "not in")
    expect_error(flagQC(obj, c("m", "nope")), "nope")
})

test_that("the absent-column error counts only what is actually absent", {
    obj <- make_metric_seurat(m = metric_values())

    expect_error(flagQC(obj, c("m", "nope")), "1 column")
    expect_error(flagQC(obj, c("nope", "also_nope")), "2 columns")
})

test_that("an invalid bound is rejected and the valid set named", {
    obj <- make_metric_seurat(m = metric_values())

    expect_error(flagQC(obj, c(m = "up")), "Invalid bound")
    expect_error(flagQC(obj, c(m = "up")), "lower")
    expect_error(flagQC(obj, c(m = "up")), "higher")
})

test_that("the same metric and bound twice is rejected", {
    # both entries write the same column, so only the second would survive
    obj <- make_metric_seurat(m = metric_values())

    expect_error(flagQC(obj, c(m = "both", m = "both")),
                 "repeats the same metric and bound")
})

test_that("the same metric with different bounds is allowed", {
    obj <- make_metric_seurat(m = metric_values(100))

    res <- flagQC(obj, c(m = "lower", m = "higher"), log = FALSE)

    expect_true("m_3_lower_outliers" %in% colnames(res@meta.data))
    expect_true("m_3_higher_outliers" %in% colnames(res@meta.data))
})

test_that("a non-numeric metric column is rejected and its class named", {
    obj <- make_metric_seurat(m = metric_values())
    obj$chr <- rep_len(letters[1:2], ncol(obj))

    expect_error(flagQC(obj, "chr"), "is not numeric")
    expect_error(flagQC(obj, "chr"), "<character>")
})

test_that("negative values under log are an error, not a silent NaN", {
    obj <- make_metric_seurat(m = metric_values() - 100)

    expect_error(flagQC(obj, "m"), "cannot be log-transformed")
    expect_error(flagQC(obj, "m"), "log = FALSE")
})

test_that("negative values are fine with log = FALSE", {
    obj <- make_metric_seurat(m = metric_values(100) - 100)
    expect_no_error(flagQC(obj, "m", log = FALSE))
})

# ---- how metrics and bounds pair up --------------------------------------

test_that("unnamed entries inherit the function-level bound", {
    obj <- make_metric_seurat(m = metric_values(100))

    res <- flagQC(obj, "m", bound = "higher", log = FALSE)

    expect_true("m_3_higher_outliers" %in% colnames(res@meta.data))
})

test_that("named entries override the function-level bound", {
    obj <- make_metric_seurat(m = metric_values(100))

    res <- flagQC(obj, c(m = "lower"), bound = "higher", log = FALSE)

    expect_true("m_3_lower_outliers" %in% colnames(res@meta.data))
    expect_false("m_3_higher_outliers" %in% colnames(res@meta.data))
})

test_that("named and unnamed entries can be mixed in one call", {
    obj <- make_metric_seurat(a = metric_values(100), b = metric_values(100))

    res <- flagQC(obj, c("a", b = "higher"), bound = "lower", log = FALSE)

    expect_true("a_3_lower_outliers" %in% colnames(res@meta.data))
    expect_true("b_3_higher_outliers" %in% colnames(res@meta.data))
})

test_that("bound itself is matched, not taken verbatim", {
    obj <- make_metric_seurat(m = metric_values(100))

    expect_error(flagQC(obj, "m", bound = "up"), "'arg' should be one of")
})

# ---- group_by ------------------------------------------------------------

test_that("an absent group_by column is rejected", {
    obj <- make_metric_seurat(m = metric_values())

    expect_error(flagQC(obj, "m", group_by = "nope"), "is not in")
})

test_that("group_by must be a single column name", {
    obj <- make_metric_seurat(m = metric_values())

    expect_error(flagQC(obj, "m", group_by = c("Ident", "m")),
                 "single metadata column name")
    expect_error(flagQC(obj, "m", group_by = 1),
                 "single metadata column name")
    expect_error(flagQC(obj, "m", group_by = NA_character_),
                 "single metadata column name")
})

test_that("a continuous group_by column is rejected", {
    # one group per cell means no thresholds at all
    obj <- make_metric_seurat(m = metric_values())

    expect_error(flagQC(obj, "m", group_by = "m"), "is continuous")
    expect_error(flagQC(obj, "m", group_by = "m"), "factor")
})

test_that("an entirely missing group_by column is rejected", {
    obj <- make_metric_seurat(m = metric_values())
    obj$grp <- rep(NA_character_, ncol(obj))

    expect_error(flagQC(obj, "m", group_by = "grp"), "entirely missing")
})

test_that("validation runs before any column is written", {
    obj <- make_metric_seurat(a = metric_values(100), b = metric_values(100))

    expect_error(flagQC(obj, c("a", "nope")), "not in")
    # the object the caller still holds is untouched either way, but the
    # failure must come from validation, not from a half-finished loop
    expect_false("a_3_both_outliers" %in% colnames(obj@meta.data))
})
