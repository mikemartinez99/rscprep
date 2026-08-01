# Single-path mode. The defining contract here (vs. manifest mode) is that a
# missing H5 file is FATAL -- there is no other sample to fall back on.

# ---- non-multiplexed -----------------------------------------------------

test_that("non-multiplexed path reads the run's H5 and names by parent dir", {
    root  <- local_root()
    calls <- local_mock_read10x()
    outs  <- make_count_run(root, "Sample_1")

    res <- seuratFrom10x(path = outs, multiplexed = FALSE)

    expect_type(res, "list")
    expect_length(res, 1)
    expect_named(res, "Sample_1")
    expect_s4_class(res$Sample_1, "Seurat")

    # the sample name comes from the directory ABOVE outs, not "outs"
    expect_equal(calls$files, file.path(outs, "filtered_feature_bc_matrix.h5"))
})

test_that("Ident, project and assay are set on the returned object", {
    root <- local_root()
    local_mock_read10x()
    outs <- make_count_run(root, "Sample_1")

    obj <- seuratFrom10x(path = outs, multiplexed = FALSE)$Sample_1

    expect_true(all(obj$Ident == "Sample_1"))
    expect_equal(as.character(unique(obj$orig.ident)), "Sample_1")
    expect_equal(Seurat::DefaultAssay(obj), "RNA")
})

test_that("trailing slashes on path do not corrupt the sample name", {
    root <- local_root()
    outs <- make_count_run(root, "Sample_1")

    for (p in c(outs, paste0(outs, "/"), paste0(outs, "///"))) {
        calls <- local_mock_read10x()
        res   <- seuratFrom10x(path = p, multiplexed = FALSE)

        expect_named(res, "Sample_1", info = p)
        # and the H5 path is still well-formed (no doubled separators)
        expect_equal(calls$files,
                     file.path(outs, "filtered_feature_bc_matrix.h5"),
                     info = p)
    }
})

test_that("missing H5 in non-multiplexed path mode is an error, not a warning", {
    root <- local_root()
    local_mock_read10x()
    outs <- make_count_run(root, "Sample_1", h5 = FALSE)

    expect_error(seuratFrom10x(path = outs, multiplexed = FALSE),
                 "File not found at")
    # the message names the exact path that was probed
    expect_error(seuratFrom10x(path = outs, multiplexed = FALSE),
                 file.path(outs, "filtered_feature_bc_matrix.h5"),
                 fixed = TRUE)
    expect_no_warning(
        try(seuratFrom10x(path = outs, multiplexed = FALSE), silent = TRUE)
    )
})

test_that("a nonexistent path errors rather than returning an empty list", {
    local_mock_read10x()
    expect_error(
        seuratFrom10x(path = "/no/such/Sample_1/outs", multiplexed = FALSE),
        "File not found at"
    )
})

# ---- multiplexed ---------------------------------------------------------

test_that("multiplexed path expands into one object per demuxed sample", {
    root  <- local_root()
    calls <- local_mock_read10x()
    outs  <- make_multi_run(root, "Pool1", c("Sample_A", "Sample_B"))

    res <- seuratFrom10x(path = outs, multiplexed = TRUE)

    expect_length(res, 2)
    expect_setequal(names(res), c("Sample_A", "Sample_B"))
    expect_true(all(vapply(res, inherits, logical(1), "Seurat")))

    # reads the per-sample filename, not the top-level one
    expect_setequal(
        calls$files,
        file.path(outs, "per_sample_outs", c("Sample_A", "Sample_B"),
                  "sample_filtered_feature_bc_matrix.h5")
    )
})

test_that("each demuxed sample carries its own Ident", {
    root <- local_root()
    local_mock_read10x()
    outs <- make_multi_run(root, "Pool1", c("Sample_A", "Sample_B"))

    res <- seuratFrom10x(path = outs, multiplexed = TRUE)

    expect_true(all(res$Sample_A$Ident == "Sample_A"))
    expect_true(all(res$Sample_B$Ident == "Sample_B"))
})

test_that("multiplexed mode does not descend into nested subdirectories", {
    root <- local_root()
    local_mock_read10x()
    outs <- make_multi_run(root, "Pool1", "Sample_A")
    dir.create(file.path(outs, "per_sample_outs", "Sample_A", "count"),
               recursive = TRUE)

    res <- seuratFrom10x(path = outs, multiplexed = TRUE)

    expect_named(res, "Sample_A")   # "count" must not appear as a sample
})

test_that("missing H5 for a demuxed sample is fatal in path mode", {
    root <- local_root()
    local_mock_read10x()
    outs <- make_multi_run(root, "Pool1", c("Sample_A", "Sample_B"),
                           h5 = c(TRUE, FALSE))

    expect_error(seuratFrom10x(path = outs, multiplexed = TRUE),
                 "File not found at")
})

test_that("multiplexed run with no per_sample_outs yields an empty list", {
    # Documents current behaviour: list.dirs() on a nonexistent directory
    # returns character(0), so this is silently empty rather than an error.
    root <- local_root()
    local_mock_read10x()
    outs <- make_multi_run(root, "Pool1", make_per_sample_outs = FALSE)

    res <- seuratFrom10x(path = outs, multiplexed = TRUE)

    expect_type(res, "list")
    expect_length(res, 0)
})

test_that("wrong 'multiplexed' flag surfaces as a missing-file error", {
    # A user who flags a plain count run as multiplexed gets an empty list;
    # the reverse gets a clear file-not-found naming the file we looked for.
    root <- local_root()
    local_mock_read10x()
    outs <- make_multi_run(root, "Pool1", "Sample_A")

    expect_error(seuratFrom10x(path = outs, multiplexed = FALSE),
                 "filtered_feature_bc_matrix\\.h5")
})

test_that("multiplexed = NA is not caught by validation", {
    # Known gap: is.null(NA) is FALSE, so the "you must specify multiplexed"
    # guard passes and the failure surfaces from `if (multiplexed)` instead.
    # Update this test if NA handling is tightened.
    root <- local_root()
    local_mock_read10x()
    outs <- make_count_run(root, "Sample_1")

    expect_error(seuratFrom10x(path = outs, multiplexed = NA),
                 "missing value where TRUE/FALSE needed")
})

# ---- filtering pass-through ---------------------------------------------

test_that("min.cells and min.features reach CreateSeuratObject", {
    root <- local_root()
    local_mock_read10x()
    outs <- make_count_run(root, "Sample_1")

    # fixture: 61 features x 21 cells, of which RareGene is in 2 cells and
    # SparseCell has 5 features
    unfiltered <- seuratFrom10x(path = outs, multiplexed = FALSE,
                                min.cells = 0, min.features = 0)$Sample_1
    expect_equal(nrow(unfiltered), 61)
    expect_equal(ncol(unfiltered), 21)
    expect_true("RareGene" %in% rownames(unfiltered))
    expect_true("SparseCell" %in% colnames(unfiltered))

    filtered <- seuratFrom10x(path = outs, multiplexed = FALSE,
                              min.cells = 3, min.features = 6)$Sample_1
    expect_false("RareGene" %in% rownames(filtered))
    expect_false("SparseCell" %in% colnames(filtered))
})

test_that("package defaults filter the marginal fixture entries", {
    root <- local_root()
    local_mock_read10x()
    outs <- make_count_run(root, "Sample_1")

    obj <- seuratFrom10x(path = outs, multiplexed = FALSE)$Sample_1

    expect_equal(nrow(obj), 60)   # RareGene dropped by min.cells = 10
    expect_equal(ncol(obj), 20)   # SparseCell dropped by min.features = 50
})
