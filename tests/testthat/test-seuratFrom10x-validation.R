# Argument validation. None of these should touch the filesystem or Seurat --
# they must fail before any I/O is attempted.

test_that("requires exactly one of 'manifest' or 'path'", {
    expect_error(seuratFrom10x(), "Provide either 'manifest' or 'path'")
    expect_error(
        seuratFrom10x(manifest = data.frame(Path = "x", multiplexed = TRUE),
                      path     = "x"),
        "Provide only one of 'manifest' or 'path'"
    )
})

test_that("'path' mode requires an explicit 'multiplexed'", {
    expect_error(seuratFrom10x(path = "/some/outs"),
                 "you must also specify 'multiplexed'")
})

test_that("'multiplexed' is not required -- and is ignored -- in manifest mode", {
    root <- local_root()
    calls <- local_mock_read10x()
    outs <- make_count_run(root, "Sample_1")

    manifest <- data.frame(Path = outs, multiplexed = FALSE)

    # `multiplexed = TRUE` here contradicts the manifest row; the row wins,
    # so the non-multiplexed H5 is the one that gets read.
    res <- seuratFrom10x(manifest = manifest, multiplexed = TRUE)

    expect_named(res, "Sample_1")
    expect_equal(calls$files, file.path(outs, "filtered_feature_bc_matrix.h5"))
})

test_that("manifest must carry the required columns", {
    expect_error(seuratFrom10x(manifest = data.frame(multiplexed = TRUE)),
                 "missing required column\\(s\\): Path")
    expect_error(seuratFrom10x(manifest = data.frame(Path = "x")),
                 "missing required column\\(s\\): multiplexed")

    # both missing -> both named, in the declared order
    expect_error(seuratFrom10x(manifest = data.frame(other = 1)),
                 "missing required column\\(s\\): Path, multiplexed")
})

test_that("manifest column check precedes the empty-manifest check", {
    # A zero-row frame that is *also* missing columns should report the
    # columns, since that is the more actionable error.
    expect_error(seuratFrom10x(manifest = data.frame(other = character(0))),
                 "missing required column\\(s\\)")
})

test_that("empty manifest is rejected", {
    empty <- data.frame(Path = character(0), multiplexed = logical(0))
    expect_error(seuratFrom10x(manifest = empty), "manifest has no rows")
})

test_that("manifest$multiplexed must be logical, not coercible-to-logical", {
    expect_error(
        seuratFrom10x(manifest = data.frame(Path = "x", multiplexed = "TRUE")),
        "must be logical"
    )
    expect_error(
        seuratFrom10x(manifest = data.frame(Path = "x", multiplexed = 1)),
        "must be logical"
    )
    expect_error(
        seuratFrom10x(manifest = data.frame(Path = "x",
                                            multiplexed = factor("TRUE"))),
        "must be logical"
    )
})

test_that("extra manifest columns are ignored", {
    root <- local_root()
    local_mock_read10x()
    outs <- make_count_run(root, "Sample_1")

    # matches the shape used in R/testing.R, which carries a Sample column
    manifest <- data.frame(Sample      = "Pool_1",
                           Path        = outs,
                           multiplexed = FALSE,
                           note        = "anything")

    expect_named(seuratFrom10x(manifest = manifest), "Sample_1")
})

test_that("a tibble-like manifest is accepted", {
    skip_if_not_installed("tibble")

    root <- local_root()
    local_mock_read10x()
    outs <- make_count_run(root, "Sample_1")

    manifest <- tibble::tibble(Path = outs, multiplexed = FALSE)
    expect_named(seuratFrom10x(manifest = manifest), "Sample_1")
})
