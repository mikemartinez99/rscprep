# Same fixture strategy as the filtered-matrix tests: real directory layouts
# with empty placeholder H5 files, and Seurat::Read10X_h5() mocked, so the
# thing under test is path resolution and naming rather than HDF5 parsing.
# What is asserted throughout is *which* file the package asked for, not
# merely that it asked for one.

test_that("path mode reads the flat raw matrix and names it after the run", {
    root  <- local_root()
    calls <- local_mock_read10x()
    outs  <- make_count_run(root, "Sample_1")
    add_raw_h5(outs, "flat")

    res <- rawMatrixFrom10x(path = outs, multi = FALSE)

    expect_named(res, "Sample_1")
    expect_equal(calls$files, file.path(outs, "raw_feature_bc_matrix.h5"))
})

test_that("a multiplexed run yields one run-level matrix, not one per sample", {
    root  <- local_root()
    calls <- local_mock_read10x()
    # two demuxed samples in per_sample_outs; the raw matrix is shared
    outs  <- make_multi_run(root, "Pool1", c("Sample_A", "Sample_B"))
    add_raw_h5(outs, "multi")

    res <- rawMatrixFrom10x(path = outs, multi = TRUE)

    expect_named(res, "Pool1")
    expect_length(res, 1)
    expect_equal(calls$files,
                 file.path(outs, "multi", "count", "raw_feature_bc_matrix.h5"))
})

test_that("a multiplexed run falls back to the flat location", {
    root  <- local_root()
    calls <- local_mock_read10x()
    outs  <- make_multi_run(root, "Pool1")
    add_raw_h5(outs, "flat")   # no multi/count/ -- e.g. staged output

    expect_named(rawMatrixFrom10x(path = outs, multi = TRUE), "Pool1")
    expect_equal(calls$files, file.path(outs, "raw_feature_bc_matrix.h5"))
})

test_that("the run-level location wins when both exist", {
    root  <- local_root()
    calls <- local_mock_read10x()
    outs  <- make_multi_run(root, "Pool1")
    add_raw_h5(outs, "flat")
    add_raw_h5(outs, "multi")

    rawMatrixFrom10x(path = outs, multi = TRUE)

    expect_equal(calls$files,
                 file.path(outs, "multi", "count", "raw_feature_bc_matrix.h5"))
})

test_that("trailing slashes on the path do not change the sample name", {
    root <- local_root()
    local_mock_read10x()
    outs <- make_count_run(root, "Sample_1")
    add_raw_h5(outs, "flat")

    expect_named(rawMatrixFrom10x(path = paste0(outs, "///"), multi = FALSE),
                 "Sample_1")
})

test_that("a missing raw matrix is fatal in path mode", {
    root <- local_root()
    local_mock_read10x()
    outs <- make_count_run(root, "Sample_1")   # filtered H5 only, no raw

    expect_error(rawMatrixFrom10x(path = outs, multi = FALSE),
                 "Raw matrix for")
})

test_that("the error names every location that was searched", {
    root <- local_root()
    local_mock_read10x()
    outs <- make_multi_run(root, "Pool1")

    expect_error(rawMatrixFrom10x(path = outs, multi = TRUE),
                 "multi/count/raw_feature_bc_matrix\\.h5")
})

test_that("manifest mode skips a missing raw matrix with a warning", {
    root <- local_root()
    local_mock_read10x()
    present <- make_count_run(root, "Sample_1")
    absent  <- make_count_run(root, "Sample_2")
    add_raw_h5(present, "flat")

    manifest <- data.frame(Path  = c(absent, present),
                           multi = c(FALSE, FALSE))

    expect_warning(res <- rawMatrixFrom10x(manifest = manifest), "Raw matrix for")
    expect_named(res, "Sample_1")
})

test_that("manifest mode mixes multiplexed and non-multiplexed rows", {
    root  <- local_root()
    calls <- local_mock_read10x()
    single <- make_count_run(root, "Sample_1")
    pool   <- make_multi_run(root, "Pool1")
    add_raw_h5(single, "flat")
    add_raw_h5(pool, "multi")

    manifest <- data.frame(Path  = c(single, pool),
                           multi = c(FALSE, TRUE))

    expect_named(rawMatrixFrom10x(manifest = manifest), c("Sample_1", "Pool1"))
    expect_equal(
        calls$files,
        c(file.path(single, "raw_feature_bc_matrix.h5"),
          file.path(pool, "multi", "count", "raw_feature_bc_matrix.h5"))
    )
})

test_that("duplicate run names are made unique, with a warning", {
    root <- local_root()
    local_mock_read10x()
    # same run directory listed twice -- e.g. a re-sequenced lane
    outs <- make_count_run(root, "Sample_1")
    add_raw_h5(outs, "flat")

    manifest <- data.frame(Path = c(outs, outs), multi = c(FALSE, FALSE))

    expect_warning(res <- rawMatrixFrom10x(manifest = manifest), "Duplicate run name")
    expect_named(res, c("Sample_1", "Sample_1_1"))
})

test_that("a multi-modality H5 is reduced to Gene Expression", {
    root <- local_root()
    local_mock_read10x(counts = list(
        "Gene Expression"        = make_counts(n_genes = 30, extras = FALSE),
        "Antibody Capture"       = make_counts(n_genes = 4,  extras = FALSE),
        "Multiplexing Capture"   = make_counts(n_genes = 2,  extras = FALSE)
    ))
    outs <- make_count_run(root, "Sample_1")
    add_raw_h5(outs, "flat")

    res <- rawMatrixFrom10x(path = outs, multi = FALSE)

    expect_equal(nrow(res$Sample_1), 30)
})

test_that("a modality list with no gene expression warns and takes the first", {
    root <- local_root()
    local_mock_read10x(counts = list(
        "Antibody Capture" = make_counts(n_genes = 4, extras = FALSE)
    ))
    outs <- make_count_run(root, "Sample_1")
    add_raw_h5(outs, "flat")

    expect_warning(res <- rawMatrixFrom10x(path = outs, multi = FALSE),
                   "no .*Gene Expression.* feature type")
    expect_equal(nrow(res$Sample_1), 4)
})

test_that("the result is the matrix itself, not a Seurat object", {
    root <- local_root()
    local_mock_read10x()
    outs <- make_count_run(root, "Sample_1")
    add_raw_h5(outs, "flat")

    res <- rawMatrixFrom10x(path = outs, multi = FALSE)

    expect_false(inherits(res$Sample_1, "Seurat"))
    expect_true(inherits(res$Sample_1, "Matrix"))
})

test_that("argument validation matches seuratFrom10x()", {
    expect_error(rawMatrixFrom10x(), "Provide either `manifest` or `path`")
    expect_error(
        rawMatrixFrom10x(manifest = data.frame(Path = "x", multi = TRUE),
                         path     = "x"),
        "Provide only one of `manifest` or `path`"
    )
    expect_error(rawMatrixFrom10x(path = "/some/outs"),
                 "you must also specify `multi`")
    expect_error(rawMatrixFrom10x(manifest = data.frame(other = 1)),
                 "missing required columns: Path and multi")
    expect_error(
        rawMatrixFrom10x(manifest = data.frame(Path = character(0),
                                               multi = logical(0))),
        "`manifest` has no rows"
    )
    expect_error(
        rawMatrixFrom10x(manifest = data.frame(Path = "x", multi = "TRUE")),
        "must be logical"
    )
})
