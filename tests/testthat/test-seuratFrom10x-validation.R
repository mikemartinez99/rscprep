# Argument validation. None of these should touch the filesystem or Seurat --
# they must fail before any I/O is attempted.

test_that("requires exactly one of 'manifest' or 'path'", {
    expect_error(seuratFrom10x(), "Provide either `manifest` or `path`")
    expect_error(
        seuratFrom10x(manifest = data.frame(Path = "x", multi = TRUE),
                      path     = "x"),
        "Provide only one of `manifest` or `path`"
    )
})

test_that("'path' mode requires an explicit 'multi'", {
    expect_error(seuratFrom10x(path = "/some/outs"),
                 "you must also specify `multi`")
})

test_that("'multi' is not required -- and is ignored -- in manifest mode", {
    root <- local_root()
    calls <- local_mock_read10x()
    outs <- make_count_run(root, "Sample_1")

    manifest <- data.frame(Path = outs, multi = FALSE)

    # `multi = TRUE` here contradicts the manifest row; the row wins,
    # so the `cellranger count` H5 is the one that gets read.
    res <- seuratFrom10x(manifest = manifest, multi = TRUE)

    expect_named(res, "Sample_1")
    expect_equal(calls$files, file.path(outs, "filtered_feature_bc_matrix.h5"))
})

test_that("manifest must carry the required columns", {
    # cli pluralises the noun, so a single missing column reads "column"
    expect_error(seuratFrom10x(manifest = data.frame(multi = TRUE)),
                 "missing required column: Path")
    expect_error(seuratFrom10x(manifest = data.frame(Path = "x")),
                 "missing required column: multi")

    # both missing -> both named, in the declared order
    expect_error(seuratFrom10x(manifest = data.frame(other = 1)),
                 "missing required columns: Path and multi")
})

test_that("manifest column check precedes the empty-manifest check", {
    # A zero-row frame that is *also* missing columns should report the
    # columns, since that is the more actionable error.
    expect_error(seuratFrom10x(manifest = data.frame(other = character(0))),
                 "missing required column")
})

test_that("empty manifest is rejected", {
    empty <- data.frame(Path = character(0), multi = logical(0))
    expect_error(seuratFrom10x(manifest = empty), "`manifest` has no rows")
})

test_that("manifest$multi must be logical, not coercible-to-logical", {
    expect_error(
        seuratFrom10x(manifest = data.frame(Path = "x", multi = "TRUE")),
        "must be logical"
    )
    expect_error(
        seuratFrom10x(manifest = data.frame(Path = "x", multi = 1)),
        "must be logical"
    )
    expect_error(
        seuratFrom10x(manifest = data.frame(Path = "x",
                                            multi = factor("TRUE"))),
        "must be logical"
    )
})

test_that("extra manifest columns are ignored", {
    root <- local_root()
    local_mock_read10x()
    outs <- make_count_run(root, "Sample_1")

    # matches the shape used in R/testing.R, which carries a Sample column
    manifest <- data.frame(Sample = "Pool_1",
                           Path   = outs,
                           multi  = FALSE,
                           note   = "anything")

    expect_named(seuratFrom10x(manifest = manifest), "Sample_1")
})

test_that("a tibble-like manifest is accepted", {
    skip_if_not_installed("tibble")

    root <- local_root()
    local_mock_read10x()
    outs <- make_count_run(root, "Sample_1")

    manifest <- tibble::tibble(Path = outs, multi = FALSE)
    expect_named(seuratFrom10x(manifest = manifest), "Sample_1")
})

# ---- deprecated `multiplexed` spelling ------------------------------------
#
# `multiplexed` was renamed to `multi` because it selects a Cell Ranger output
# layout rather than describing the experiment. The old spelling still works,
# so manifests and scripts written against the previous name keep running.

test_that("the `multiplexed` argument still works, with a warning", {
    root  <- local_root()
    calls <- local_mock_read10x()
    outs  <- make_count_run(root, "Sample_1")

    expect_warning(res <- seuratFrom10x(path = outs, multiplexed = FALSE),
                   "deprecated")

    expect_named(res, "Sample_1")
    expect_equal(calls$files, file.path(outs, "filtered_feature_bc_matrix.h5"))
})

test_that("the deprecated argument still selects the multi layout", {
    root  <- local_root()
    calls <- local_mock_read10x()
    outs  <- make_multi_run(root, "Pool1", "Sample_A")

    expect_warning(res <- seuratFrom10x(path = outs, multiplexed = TRUE),
                   "deprecated")

    expect_named(res, "Sample_A")
    expect_equal(calls$files,
                 file.path(outs, "per_sample_outs", "Sample_A", "count",
                           "sample_filtered_feature_bc_matrix.h5"))
})

test_that("supplying both spellings is an error", {
    expect_error(
        seuratFrom10x(path = "/some/outs", multi = TRUE, multiplexed = TRUE),
        "Provide only one of `multi` or `multiplexed`"
    )
})

test_that("a manifest still using the `multiplexed` column works, with a warning", {
    root  <- local_root()
    calls <- local_mock_read10x()
    outs  <- make_count_run(root, "Sample_1")

    expect_warning(
        res <- seuratFrom10x(manifest = data.frame(Path        = outs,
                                                   multiplexed = FALSE)),
        "deprecated"
    )

    expect_named(res, "Sample_1")
    expect_equal(calls$files, file.path(outs, "filtered_feature_bc_matrix.h5"))
})

test_that("a `multi` column wins over a stale `multiplexed` column", {
    # A half-renamed manifest must not silently resolve to the old column.
    root  <- local_root()
    calls <- local_mock_read10x()
    outs  <- make_multi_run(root, "Pool1", "Sample_A")

    manifest <- data.frame(Path        = outs,
                           multi       = TRUE,
                           multiplexed = FALSE)

    expect_named(seuratFrom10x(manifest = manifest), "Sample_A")
    expect_equal(calls$files,
                 file.path(outs, "per_sample_outs", "Sample_A", "count",
                           "sample_filtered_feature_bc_matrix.h5"))
})

test_that("the deprecated column is validated like the new one", {
    expect_warning(
        expect_error(
            seuratFrom10x(manifest = data.frame(Path = "x", multiplexed = "TRUE")),
            "must be logical"
        ),
        "deprecated"
    )
})
