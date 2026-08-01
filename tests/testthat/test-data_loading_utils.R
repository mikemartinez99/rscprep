# Internal helpers. Tests run inside the package namespace, so the @noRd
# functions are reachable directly.
#
# The behaviour worth pinning down here is the `on_missing` callback: it is
# injected rather than hardcoded, and it is the single mechanism implementing
# "fatal in path mode, skip-with-warning in manifest mode".

# ---- readOneSample10x ----------------------------------------------------

test_that("readOneSample10x builds an object from the expected file", {
    root  <- local_root()
    calls <- local_mock_read10x()
    outs  <- make_count_run(root, "Sample_1")

    obj <- readOneSample10x(outs, "Sample_1", min.cells = 0, min.features = 0)

    expect_s4_class(obj, "Seurat")
    expect_equal(calls$files, file.path(outs, "filtered_feature_bc_matrix.h5"))
    expect_true(all(obj$Ident == "Sample_1"))
    expect_equal(as.character(unique(obj$orig.ident)), "Sample_1")
})

test_that("readOneSample10x honours the filename argument", {
    root  <- local_root()
    calls <- local_mock_read10x()
    outs  <- make_multi_run(root, "Pool1", "Sample_A")
    sdir  <- file.path(outs, "per_sample_outs", "Sample_A")

    readOneSample10x(sdir, "Sample_A", 0, 0,
                     filename = "sample_filtered_feature_bc_matrix.h5")

    expect_equal(calls$files,
                 file.path(sdir, "sample_filtered_feature_bc_matrix.h5"))
})

test_that("Ident is independent of the directory name", {
    # readOneSample10x takes ident as an argument; it must not re-derive it
    root <- local_root()
    local_mock_read10x()
    outs <- make_count_run(root, "Sample_1")

    obj <- readOneSample10x(outs, "CustomLabel", 0, 0)

    expect_true(all(obj$Ident == "CustomLabel"))
})

test_that("readOneSample10x defaults to warning and returns NULL when missing", {
    root <- local_root()
    calls <- local_mock_read10x()
    outs <- make_count_run(root, "Sample_1", h5 = FALSE)

    expect_warning(obj <- readOneSample10x(outs, "Sample_1", 0, 0),
                   "File not found at")
    expect_null(obj)
    # and it short-circuits before attempting to read
    expect_length(calls$files, 0)
})

test_that("readOneSample10x routes the missing-file message through on_missing", {
    root <- local_root()
    local_mock_read10x()
    outs <- make_count_run(root, "Sample_1", h5 = FALSE)

    # fatal variant, as seuratFrom10x() supplies in path mode
    expect_error(
        readOneSample10x(outs, "Sample_1", 0, 0,
                         on_missing = function(m) stop(m, call. = FALSE)),
        "File not found at"
    )

    # arbitrary handler: message is passed as a single string, and control
    # returns to the caller (NULL) if the handler does not abort
    seen <- NULL
    obj <- readOneSample10x(outs, "Sample_1", 0, 0,
                            on_missing = function(m) seen <<- m)

    expect_type(seen, "character")
    expect_length(seen, 1)
    expect_match(seen, file.path(outs, "filtered_feature_bc_matrix.h5"),
                 fixed = TRUE)
    expect_null(obj)
})

# ---- readMultiplexedRun10x -----------------------------------------------

test_that("readMultiplexedRun10x names entries after per_sample_outs dirs", {
    root <- local_root()
    local_mock_read10x()
    outs <- make_multi_run(root, "Pool1", c("Sample_A", "Sample_B"))

    res <- readMultiplexedRun10x(outs, 0, 0)

    expect_setequal(names(res), c("Sample_A", "Sample_B"))
    expect_true(all(vapply(res, inherits, logical(1), "Seurat")))
    expect_true(all(res$Sample_A$Ident == "Sample_A"))
})

test_that("readMultiplexedRun10x reads the per-sample filename", {
    root  <- local_root()
    calls <- local_mock_read10x()
    outs  <- make_multi_run(root, "Pool1", "Sample_A")

    readMultiplexedRun10x(outs, 0, 0)

    expect_equal(calls$files,
                 file.path(outs, "per_sample_outs", "Sample_A",
                           "sample_filtered_feature_bc_matrix.h5"))
})

test_that("readMultiplexedRun10x keeps a NULL slot for a missing sample", {
    # The helper itself does not drop failures -- it returns a same-length,
    # fully-named list with NULL in the failed slot. seuratFrom10x() is what
    # filters NULLs out after unlisting.
    root <- local_root()
    local_mock_read10x()
    outs <- make_multi_run(root, "Pool1", c("Sample_A", "Sample_B"),
                           h5 = c(TRUE, FALSE))

    expect_warning(res <- readMultiplexedRun10x(outs, 0, 0),
                   "File not found at")

    expect_length(res, 2)
    expect_setequal(names(res), c("Sample_A", "Sample_B"))
    expect_s4_class(res$Sample_A, "Seurat")
    expect_null(res$Sample_B)
})

test_that("readMultiplexedRun10x propagates on_missing to each sample", {
    root <- local_root()
    local_mock_read10x()
    outs <- make_multi_run(root, "Pool1", c("Sample_A", "Sample_B"),
                           h5 = c(TRUE, FALSE))

    expect_error(
        readMultiplexedRun10x(outs, 0, 0,
                              on_missing = function(m) stop(m, call. = FALSE)),
        "File not found at"
    )
})

test_that("readMultiplexedRun10x returns an empty named list when there are no samples", {
    root <- local_root()
    local_mock_read10x()
    outs <- make_multi_run(root, "Pool1", make_per_sample_outs = FALSE)

    res <- readMultiplexedRun10x(outs, 0, 0)

    expect_type(res, "list")
    expect_length(res, 0)
})

test_that("readMultiplexedRun10x ignores files sitting in per_sample_outs", {
    root <- local_root()
    local_mock_read10x()
    outs <- make_multi_run(root, "Pool1", "Sample_A")
    file.create(file.path(outs, "per_sample_outs", "metrics_summary.csv"))

    res <- readMultiplexedRun10x(outs, 0, 0)

    expect_named(res, "Sample_A")
})

test_that("readMultiplexedRun10x forwards min.cells and min.features", {
    root <- local_root()
    local_mock_read10x()
    outs <- make_multi_run(root, "Pool1", "Sample_A")

    unfiltered <- readMultiplexedRun10x(outs, 0, 0)$Sample_A
    filtered   <- readMultiplexedRun10x(outs, 3, 6)$Sample_A

    expect_equal(dim(unfiltered), c(61L, 21L))
    expect_equal(dim(filtered),   c(60L, 20L))
})
