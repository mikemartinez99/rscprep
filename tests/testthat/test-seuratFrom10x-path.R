# Single-path mode. The defining contract here (vs. manifest mode) is that a
# missing H5 file is FATAL -- there is no other sample to fall back on.

# ---- cellranger count ----------------------------------------------------

test_that("a count run reads the run's H5 and names by parent dir", {
    root  <- local_root()
    calls <- local_mock_read10x()
    outs  <- make_count_run(root, "Sample_1")

    res <- seuratFrom10x(path = outs, multi = FALSE)

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

    obj <- seuratFrom10x(path = outs, multi = FALSE)$Sample_1

    expect_true(all(obj$Ident == "Sample_1"))
    expect_equal(as.character(unique(obj$orig.ident)), "Sample_1")
    expect_equal(Seurat::DefaultAssay(obj), "RNA")
})

test_that("trailing slashes on path do not corrupt the sample name", {
    root <- local_root()
    outs <- make_count_run(root, "Sample_1")

    for (p in c(outs, paste0(outs, "/"), paste0(outs, "///"))) {
        calls <- local_mock_read10x()
        res   <- seuratFrom10x(path = p, multi = FALSE)

        expect_named(res, "Sample_1", info = p)
        # and the H5 path is still well-formed (no doubled separators)
        expect_equal(calls$files,
                     file.path(outs, "filtered_feature_bc_matrix.h5"),
                     info = p)
    }
})

test_that("missing H5 in count path mode is an error, not a warning", {
    root <- local_root()
    local_mock_read10x()
    outs <- make_count_run(root, "Sample_1", h5 = FALSE)

    expect_error(seuratFrom10x(path = outs, multi = FALSE),
                 "File not found at")
    # the message names the exact path that was probed
    expect_error(seuratFrom10x(path = outs, multi = FALSE),
                 file.path(outs, "filtered_feature_bc_matrix.h5"),
                 fixed = TRUE)
    expect_no_warning(
        try(seuratFrom10x(path = outs, multi = FALSE), silent = TRUE)
    )
})

test_that("a nonexistent path errors rather than returning an empty list", {
    local_mock_read10x()
    expect_error(
        seuratFrom10x(path = "/no/such/Sample_1/outs", multi = FALSE),
        "File not found at"
    )
})

# ---- cellranger multi ----------------------------------------------------

test_that("a multi path expands into one object per demuxed sample", {
    root  <- local_root()
    calls <- local_mock_read10x()
    outs  <- make_multi_run(root, "Pool1", c("Sample_A", "Sample_B"))

    res <- seuratFrom10x(path = outs, multi = TRUE)

    expect_length(res, 2)
    expect_setequal(names(res), c("Sample_A", "Sample_B"))
    expect_true(all(vapply(res, inherits, logical(1), "Seurat")))

    # reads each sample's H5 from its own count/ folder, not the top-level one
    expect_setequal(
        calls$files,
        file.path(outs, "per_sample_outs", c("Sample_A", "Sample_B"),
                  "count", "sample_filtered_feature_bc_matrix.h5")
    )
})

test_that("each demuxed sample carries its own Ident", {
    root <- local_root()
    local_mock_read10x()
    outs <- make_multi_run(root, "Pool1", c("Sample_A", "Sample_B"))

    res <- seuratFrom10x(path = outs, multi = TRUE)

    expect_true(all(res$Sample_A$Ident == "Sample_A"))
    expect_true(all(res$Sample_B$Ident == "Sample_B"))
})

test_that("only per_sample_outs children become samples", {
    # per_sample_outs/<sample>/ is the sample level; anything deeper (count/,
    # vdj_b/, analysis/) belongs to that sample and must not surface as one.
    root <- local_root()
    local_mock_read10x()
    outs <- make_multi_run(root, "Pool1", "Sample_A")
    dir.create(file.path(outs, "per_sample_outs", "Sample_A", "vdj_b"),
               recursive = TRUE)

    res <- seuratFrom10x(path = outs, multi = TRUE)

    expect_named(res, "Sample_A")
})

test_that("a flattened per-sample H5 is found as a fallback", {
    # Staging pipelines that strip the count/ folder still resolve.
    root  <- local_root()
    calls <- local_mock_read10x()
    outs  <- make_multi_run(root, "Pool1", "Sample_A", where = "flat")

    res <- seuratFrom10x(path = outs, multi = TRUE)

    expect_named(res, "Sample_A")
    expect_equal(calls$files,
                 file.path(outs, "per_sample_outs", "Sample_A",
                           "sample_filtered_feature_bc_matrix.h5"))
})

test_that("the count/ location wins when both are present", {
    root  <- local_root()
    calls <- local_mock_read10x()
    outs  <- make_multi_run(root, "Pool1", "Sample_A")
    file.create(sampleFilteredH5(
        file.path(outs, "per_sample_outs", "Sample_A"), "flat"))

    seuratFrom10x(path = outs, multi = TRUE)

    expect_equal(calls$files,
                 file.path(outs, "per_sample_outs", "Sample_A", "count",
                           "sample_filtered_feature_bc_matrix.h5"))
})

test_that("missing H5 for a demuxed sample is fatal in path mode", {
    root <- local_root()
    local_mock_read10x()
    outs <- make_multi_run(root, "Pool1", c("Sample_A", "Sample_B"),
                           h5 = c(TRUE, FALSE))

    expect_error(seuratFrom10x(path = outs, multi = TRUE),
                 "File not found at")
})

test_that("a run with no per_sample_outs is fatal in path mode", {
    # list.dirs() returns character(0) for an absent directory, which would
    # otherwise make the run vanish without a word.
    root <- local_root()
    local_mock_read10x()
    outs <- make_multi_run(root, "Pool1", make_per_sample_outs = FALSE)

    expect_error(seuratFrom10x(path = outs, multi = TRUE),
                 "No demuxed sample directories found")
})

test_that("an empty per_sample_outs is fatal in path mode", {
    root <- local_root()
    local_mock_read10x()
    outs <- make_multi_run(root, "Pool1", make_per_sample_outs = FALSE)
    dir.create(file.path(outs, "per_sample_outs"))

    expect_error(seuratFrom10x(path = outs, multi = TRUE),
                 "No demuxed sample directories found")
})

test_that("either wrong 'multi' flag names what was looked for", {
    root <- local_root()
    local_mock_read10x()

    # multi = FALSE on a `cellranger multi` run
    pool <- make_multi_run(root, "Pool1", "Sample_A")
    expect_error(seuratFrom10x(path = pool, multi = FALSE),
                 "filtered_feature_bc_matrix\\.h5")

    # multi = TRUE on a `cellranger count` run
    counted <- make_count_run(root, "Sample_1")
    expect_error(seuratFrom10x(path = counted, multi = TRUE),
                 "per_sample_outs")
})

test_that("multi = NA is not caught by validation", {
    # Known gap: is.null(NA) is FALSE, so the "you must specify multi"
    # guard passes and the failure surfaces from `if (multi)` instead.
    # Update this test if NA handling is tightened.
    root <- local_root()
    local_mock_read10x()
    outs <- make_count_run(root, "Sample_1")

    expect_error(seuratFrom10x(path = outs, multi = NA),
                 "missing value where TRUE/FALSE needed")
})

# ---- filtering pass-through ---------------------------------------------

test_that("min.cells and min.features reach CreateSeuratObject", {
    root <- local_root()
    local_mock_read10x()
    outs <- make_count_run(root, "Sample_1")

    # fixture: 61 features x 21 cells, of which RareGene is in 2 cells and
    # SparseCell has 5 features
    unfiltered <- seuratFrom10x(path = outs, multi = FALSE,
                                min.cells = 0, min.features = 0)$Sample_1
    expect_equal(nrow(unfiltered), 61)
    expect_equal(ncol(unfiltered), 21)
    expect_true("RareGene" %in% rownames(unfiltered))
    expect_true("SparseCell" %in% colnames(unfiltered))

    filtered <- seuratFrom10x(path = outs, multi = FALSE,
                              min.cells = 3, min.features = 6)$Sample_1
    expect_false("RareGene" %in% rownames(filtered))
    expect_false("SparseCell" %in% colnames(filtered))
})

test_that("package defaults filter the marginal fixture entries", {
    root <- local_root()
    local_mock_read10x()
    outs <- make_count_run(root, "Sample_1")

    obj <- seuratFrom10x(path = outs, multi = FALSE)$Sample_1

    expect_equal(nrow(obj), 60)   # RareGene dropped by min.cells = 10
    expect_equal(ncol(obj), 20)   # SparseCell dropped by min.features = 50
})
