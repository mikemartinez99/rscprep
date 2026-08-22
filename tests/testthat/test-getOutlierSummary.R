# getOutlierSummary() reads back the table flagQC() wrote. It is a thin
# accessor, so what is worth pinning down is the states where there is nothing
# sensible to return.
#
# Regexes use \\s+ where a message could wrap, since cli hard-wraps at the
# console width.

flagged_object <- function(...) {
    obj <- make_metric_seurat(m = metric_values(100), ...)
    flagQC(obj, c(m = "both"), log = FALSE)
}

# ---- reading the table ---------------------------------------------------

test_that("the stored table is returned", {
    res <- flagged_object()

    expect_identical(getOutlierSummary(res), res@misc$flagQC)
})

test_that("a data frame comes back", {
    qc <- getOutlierSummary(flagged_object())

    expect_s3_class(qc, "data.frame")
    expect_equal(nrow(qc), 1)
    expect_true("n_flagged" %in% names(qc))
})

test_that("accumulated rows all come back", {
    obj <- make_metric_seurat(a = metric_values(100), b = metric_values(100))
    res <- flagQC(obj, "a", log = FALSE)
    res <- flagQC(res, "b", log = FALSE)

    expect_equal(nrow(getOutlierSummary(res)), 2)
})

test_that("grouped rows come back with their group labels", {
    obj <- make_metric_seurat(m = c(metric_values(100), metric_values(100)),
                              ident = rep(c("A", "B"), each = 21))
    res <- flagQC(obj, "m", group_by = "Ident", log = FALSE)

    qc <- getOutlierSummary(res)

    expect_equal(qc$group, c("A", "B"))
})

test_that("the returned table is a plain copy, not a live view", {
    res <- flagged_object()

    qc <- getOutlierSummary(res)
    qc$n_flagged <- 999L

    expect_false(getOutlierSummary(res)$n_flagged == 999L)
})

# ---- nothing to return ---------------------------------------------------

test_that("an unflagged object is an error, not a NULL", {
    obj <- make_metric_seurat(m = metric_values(100))

    expect_error(getOutlierSummary(obj), "no outlier summary")
    expect_error(getOutlierSummary(obj), "flagQC")
})

test_that("the lookup is exact, not partial", {
    # `$` on a list partial-matches, so a differently named misc entry must
    # not be handed back as though it were the summary
    obj <- make_metric_seurat(m = metric_values(100))
    obj@misc[["flagQCsomethingelse"]] <- data.frame(x = 1)

    expect_error(getOutlierSummary(obj), "no outlier summary")
})

test_that("a misc entry that is not a data frame is an error", {
    res <- flagged_object()
    res@misc[["flagQC"]] <- "clobbered"

    expect_error(getOutlierSummary(res), "not a data frame")
    expect_error(getOutlierSummary(res), "<character>")
})

# ---- validation ----------------------------------------------------------

test_that("a list of objects is rejected with a usable hint", {
    obj_list <- make_seurat_list()

    expect_error(getOutlierSummary(obj_list), "is a list")
    expect_error(getOutlierSummary(obj_list),
                 "lapply\\(obj_list, getOutlierSummary")
})

test_that("non-Seurat input is rejected", {
    expect_error(getOutlierSummary(42), "must be a Seurat object")
})

test_that("errors name the argument as the caller wrote it", {
    my_object <- make_metric_seurat(m = metric_values(100))

    expect_error(getOutlierSummary(my_object), "my_object")
})
