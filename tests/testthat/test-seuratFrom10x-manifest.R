# Manifest mode. The defining contract here (vs. path mode) is that a missing
# H5 file WARNS and drops that sample -- one bad row must not kill the batch.

test_that("a single count row behaves like path mode", {
    root  <- local_root()
    calls <- local_mock_read10x()
    outs  <- make_count_run(root, "Sample_1")

    res <- seuratFrom10x(manifest = data.frame(Path  = outs,
                                               multi = FALSE))

    expect_named(res, "Sample_1")
    expect_s4_class(res$Sample_1, "Seurat")
    expect_equal(calls$files, file.path(outs, "filtered_feature_bc_matrix.h5"))
})

test_that("count and multi rows resolve independently", {
    root <- local_root()
    local_mock_read10x()

    s1   <- make_count_run(root, "Sample_1")
    s2   <- make_count_run(root, "Sample_2")
    pool <- make_multi_run(root, "Pool1", c("Sample_A", "Sample_B"))

    res <- seuratFrom10x(manifest = data.frame(
        Path  = c(s1, s2, pool),
        multi = c(FALSE, FALSE, TRUE)
    ))

    # the multi row expands to 2 entries -> 4 objects from 3 rows
    expect_length(res, 4)
    expect_setequal(names(res),
                    c("Sample_1", "Sample_2", "Sample_A", "Sample_B"))
    expect_true(all(vapply(res, inherits, logical(1), "Seurat")))
})

test_that("manifest row order is preserved in the returned list", {
    root <- local_root()
    local_mock_read10x()

    s1   <- make_count_run(root, "Sample_1")
    pool <- make_multi_run(root, "Pool1", c("Sample_A", "Sample_B"))
    s2   <- make_count_run(root, "Sample_2")

    res <- seuratFrom10x(manifest = data.frame(
        Path  = c(s1, pool, s2),
        multi = c(FALSE, TRUE, FALSE)
    ))

    expect_equal(names(res),
                 c("Sample_1", "Sample_A", "Sample_B", "Sample_2"))
})

test_that("each object keeps the Ident of its own sample", {
    root <- local_root()
    local_mock_read10x()

    s1   <- make_count_run(root, "Sample_1")
    pool <- make_multi_run(root, "Pool1", "Sample_A")

    res <- seuratFrom10x(manifest = data.frame(Path  = c(s1, pool),
                                               multi = c(FALSE, TRUE)))

    expect_true(all(res$Sample_1$Ident == "Sample_1"))
    expect_true(all(res$Sample_A$Ident == "Sample_A"))
})

# ---- missing-file tolerance ---------------------------------------------

test_that("a missing H5 warns and drops only that row", {
    root <- local_root()
    local_mock_read10x()

    ok      <- make_count_run(root, "Sample_1")
    missing <- make_count_run(root, "Sample_2", h5 = FALSE)

    expect_warning(
        res <- seuratFrom10x(manifest = data.frame(
            Path  = c(ok, missing),
            multi = c(FALSE, FALSE)
        )),
        "File not found at"
    )

    expect_named(res, "Sample_1")
    expect_length(res, 1)
})

test_that("a missing demuxed H5 drops one sample but keeps its siblings", {
    root <- local_root()
    local_mock_read10x()

    pool <- make_multi_run(root, "Pool1", c("Sample_A", "Sample_B"),
                           h5 = c(TRUE, FALSE))

    expect_warning(
        res <- seuratFrom10x(manifest = data.frame(Path  = pool,
                                                   multi = TRUE)),
        "File not found at"
    )

    expect_named(res, "Sample_A")
})

test_that("every row missing yields an empty list rather than an error", {
    root <- local_root()
    local_mock_read10x()

    m1 <- make_count_run(root, "Sample_1", h5 = FALSE)
    m2 <- make_count_run(root, "Sample_2", h5 = FALSE)

    suppressWarnings(
        res <- seuratFrom10x(manifest = data.frame(
            Path  = c(m1, m2),
            multi = c(FALSE, FALSE)
        ))
    )

    expect_type(res, "list")
    expect_length(res, 0)
})

test_that("a nonexistent manifest path warns instead of erroring", {
    root <- local_root()
    local_mock_read10x()
    ok <- make_count_run(root, "Sample_1")

    expect_warning(
        res <- seuratFrom10x(manifest = data.frame(
            Path  = c(ok, "/no/such/Sample_9/outs"),
            multi = c(FALSE, FALSE)
        )),
        "File not found at"
    )

    expect_named(res, "Sample_1")
})

test_that("a multi row with no per_sample_outs warns and contributes nothing", {
    # The batch contract: one broken run is reported and skipped, and the
    # rows around it still load.
    root <- local_root()
    local_mock_read10x()

    ok    <- make_count_run(root, "Sample_1")
    empty <- make_multi_run(root, "Pool1", make_per_sample_outs = FALSE)

    expect_warning(
        res <- seuratFrom10x(manifest = data.frame(Path  = c(ok, empty),
                                                   multi = c(FALSE, TRUE))),
        "No demuxed sample directories found"
    )

    expect_named(res, "Sample_1")
})

# ---- path handling -------------------------------------------------------

test_that("trailing slashes in manifest paths are stripped per row", {
    root <- local_root()
    calls <- local_mock_read10x()

    s1   <- make_count_run(root, "Sample_1")
    pool <- make_multi_run(root, "Pool1", "Sample_A")

    res <- seuratFrom10x(manifest = data.frame(
        Path  = c(paste0(s1, "/"), paste0(pool, "///")),
        multi = c(FALSE, TRUE)
    ))

    expect_setequal(names(res), c("Sample_1", "Sample_A"))
    expect_setequal(
        calls$files,
        c(file.path(s1, "filtered_feature_bc_matrix.h5"),
          file.path(pool, "per_sample_outs", "Sample_A", "count",
                    "sample_filtered_feature_bc_matrix.h5"))
    )
})

test_that("duplicate sample names across rows are both retained", {
    # Documents current behaviour: names are not deduplicated, so two runs
    # producing the same sample name yield two same-named list entries and
    # res$Name silently resolves to the first.
    root <- local_root()
    local_mock_read10x()

    a <- make_count_run(file.path(root, "batch1"), "Sample_1")
    b <- make_count_run(file.path(root, "batch2"), "Sample_1")

    res <- seuratFrom10x(manifest = data.frame(Path  = c(a, b),
                                               multi = c(FALSE, FALSE)))

    expect_length(res, 2)
    expect_equal(names(res), c("Sample_1", "Sample_1"))
})

# ---- pass-through --------------------------------------------------------

test_that("min.cells and min.features apply uniformly across all rows", {
    root <- local_root()
    local_mock_read10x()

    s1   <- make_count_run(root, "Sample_1")
    pool <- make_multi_run(root, "Pool1", "Sample_A")

    res <- seuratFrom10x(
        manifest     = data.frame(Path = c(s1, pool),
                                  multi = c(FALSE, TRUE)),
        min.cells    = 0,
        min.features = 0
    )

    expect_equal(unname(vapply(res, nrow, integer(1))), c(61L, 61L))
    expect_equal(unname(vapply(res, ncol, integer(1))), c(21L, 21L))
})

test_that("NA in manifest$multi is not caught by validation", {
    # Known gap: is.logical(NA) is TRUE, so the type check passes and the
    # failure surfaces from `if (manifest$multi[row_i])` instead.
    root <- local_root()
    local_mock_read10x()
    outs <- make_count_run(root, "Sample_1")

    expect_error(
        seuratFrom10x(manifest = data.frame(Path = outs, multi = NA)),
        "missing value where TRUE/FALSE needed"
    )
})
