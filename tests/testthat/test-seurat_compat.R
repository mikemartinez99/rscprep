# The Seurat 4 / 5 compatibility helpers. Version-neutral tests run under
# either; the ones about split layers only mean anything under Seurat 5,
# whose merge() is the only thing that creates them, so they skip under 4.

test_that("isSeuratV5() agrees with the installed SeuratObject", {
    expect_identical(isSeuratV5(),
                     utils::packageVersion("SeuratObject") >= "5.0.0")
})

test_that("getLayerValues() reads counts with genes as rows and cells as columns", {
    obj <- make_seurat("Sample_A", n_cells = 6, n_genes = 8)
    counts <- getLayerValues(obj, "counts")

    expect_identical(dimnames(counts), dimnames(obj))
    expect_true(all(counts == 1))
})

test_that("getLayerValues() returns NULL for a data layer Seurat 5 never made", {
    skip_if_not(isSeuratV5(), "Seurat 4 always fills the data slot")
    obj <- make_cellcycle_seurat(normalise = FALSE)

    expect_null(getLayerValues(obj, "data"))
})

test_that("getLayerValues() binds split layers back into every cell", {
    skip_if_not(isSeuratV5(), "only Seurat 5 splits layers on merge")
    objs   <- make_seurat_list(c("Sample_A", "Sample_B"), n_cells = 5)
    merged <- merge(objs[[1]], objs[-1], add.cell.id = names(objs))

    counts <- getLayerValues(merged, "counts")
    expect_equal(ncol(counts), 10)
    expect_setequal(colnames(counts), colnames(merged))
})

test_that("joinAssayLayers() leaves one counts layer after a Seurat 5 merge", {
    skip_if_not(isSeuratV5(), "only Seurat 5 splits layers on merge")
    objs   <- make_seurat_list(c("Sample_A", "Sample_B"))
    merged <- merge(objs[[1]], objs[-1], add.cell.id = names(objs))
    layers <- getExportedValue("SeuratObject", "Layers")

    expect_length(layers(merged), 2)
    expect_identical(layers(joinAssayLayers(merged)), "counts")
})

test_that("joinAssayLayers() is a no-op under Seurat 4", {
    skip_if(isSeuratV5(), "Seurat 4 behaviour")
    obj <- make_seurat("Sample_A")

    expect_identical(joinAssayLayers(obj), obj)
})

test_that("mergeSeurat() joins layers by default and keeps them split on request", {
    skip_if_not(isSeuratV5(), "only Seurat 5 splits layers on merge")
    objs   <- make_seurat_list(c("Sample_A", "Sample_B", "Sample_C"))
    layers <- getExportedValue("SeuratObject", "Layers")

    joined <- suppressMessages(mergeSeurat(objs))
    split  <- suppressMessages(mergeSeurat(objs, joinLayers = FALSE))

    expect_identical(layers(joined), "counts")
    expect_setequal(layers(split), paste0("counts.", names(objs)))
})

test_that("runCellCycle() normalises raw data under either version", {
    obj <- make_cellcycle_seurat(normalise = FALSE)

    expect_true(looksUnnormalised(obj))
    res <- suppressWarnings(runCellCycle(obj, "hsa"))
    expect_false(looksUnnormalised(res))
    expect_true(all(c("S.Score", "G2M.Score", "Phase") %in% colnames(res@meta.data)))
})
