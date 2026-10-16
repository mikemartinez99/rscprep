# The only tests that read real H5 files: everything else mocks Read10X_h5().
# These check the bundled example data still round-trips through the readers,
# and that the object in data/ matches the files it was built from.

skip_if_not_installed("hdf5r")

outs <- system.file("extdata", "heart_1k_v3", "outs", package = "rscprep")

test_that("bundled outs/ reads through seuratFrom10x()", {
    seurat_list <- suppressMessages(
        seuratFrom10x(path = outs, multi = FALSE, min.cells = 0, min.features = 0)
    )
    expect_named(seurat_list, "heart_1k_v3")
    # Compared through dimnames and per-cell totals rather than the counts
    # matrix itself, which SeuratObject 4 and 5 fetch with different arguments.
    obj <- seurat_list[[1]]
    expect_identical(dimnames(obj), dimnames(heart1k))
    expect_equal(unname(obj$nCount_RNA), unname(heart1k$nCount_RNA))
    expect_equal(unname(obj$nFeature_RNA), unname(heart1k$nFeature_RNA))
})

test_that("bundled raw matrix is a superset of the filtered cells", {
    raw <- suppressMessages(rawMatrixFrom10x(path = outs, multi = FALSE))[[1]]
    expect_true(all(colnames(heart1k) %in% colnames(raw)))
    expect_identical(rownames(raw), rownames(heart1k))
})

test_that("heart1k keeps the genes the QC functions look for", {
    expect_equal(ncol(heart1k), 1011)
    expect_length(grep("^mt", rownames(heart1k)), 13)
    expect_gt(length(grep("^Rp[sl]", rownames(heart1k))), 0)
})
