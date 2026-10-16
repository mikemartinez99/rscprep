# mergeSeurat() collapses the named list seuratFrom10x() returns into a single
# object. Most of its surface is validation, and that is deliberate: the list
# names become cell barcode prefixes, so a malformed list has to be rejected
# up front rather than merged into something quietly wrong.
#
# Regexes below use \\s+ where a message could wrap, since cli hard-wraps at
# the console width.

# ---- empty and single-element input --------------------------------------

test_that("an empty list is an error", {
    obj_list <- list()
    expect_error(mergeSeurat(obj_list), "is empty")
})

test_that("errors name the argument as the caller wrote it", {
    my_samples <- list()
    expect_error(mergeSeurat(my_samples), "my_samples")
})

test_that("a single-element list is returned unmerged", {
    obj_list <- make_seurat_list("Sample_A")

    res <- mergeSeurat(obj_list)

    expect_s4_class(res, "Seurat")
    expect_identical(res, obj_list$Sample_A)
    # nothing was merged, so barcodes are NOT prefixed
    expect_equal(colnames(res), paste0("Cell", 1:10))
})

test_that("the single-element shortcut runs before any validation", {
    # Documents current behaviour: a one-element list is handed back as-is
    # without ever being checked, so this returns the string unharmed.
    obj_list <- list(Sample_A = "not a seurat object")

    expect_identical(mergeSeurat(obj_list), "not a seurat object")
})

# ---- element validity ----------------------------------------------------

test_that("non-Seurat elements are rejected and named", {
    obj_list <- list(Sample_A = make_seurat("Sample_A"),
                     Sample_B = "not a seurat object")

    expect_error(mergeSeurat(obj_list), "not valid Seurat objects")
    expect_error(mergeSeurat(obj_list), "Sample_B")
})

test_that("objects that inherit from Seurat but fail validObject are rejected", {
    # inherits() alone would pass this; the validObject() half of the check is
    # what catches it.
    obj_list <- list(Sample_A = make_seurat("Sample_A"),
                     Sample_B = make_fake_seurat())

    expect_true(inherits(obj_list$Sample_B, "Seurat"))
    expect_error(mergeSeurat(obj_list), "not valid Seurat objects")
})

test_that("a bad element in an unnamed list is reported by position", {
    obj_list <- list(make_seurat("Sample_A"), "not a seurat object")

    expect_error(mergeSeurat(obj_list), "not valid Seurat objects:\\s+2")
})

test_that("validity is checked before names", {
    # unnamed AND invalid -- the invalid-object error is the one raised
    obj_list <- list("not a seurat object", "also not")

    expect_error(mergeSeurat(obj_list), "not valid Seurat objects")
})

# ---- names ---------------------------------------------------------------

test_that("an unnamed list is rejected", {
    obj_list <- unname(make_seurat_list(c("Sample_A", "Sample_B")))

    expect_error(mergeSeurat(obj_list), "must be fully named")
})

test_that("blank names are rejected", {
    obj_list <- make_seurat_list(c("Sample_A", "Sample_B"))
    names(obj_list) <- c("Sample_A", "")

    expect_error(mergeSeurat(obj_list), "must be fully named")
})

test_that("NA names are rejected", {
    obj_list <- make_seurat_list(c("Sample_A", "Sample_B"))
    names(obj_list) <- c("Sample_A", NA)

    expect_error(mergeSeurat(obj_list), "must be fully named")
})

test_that("duplicate names are rejected and listed", {
    obj_list <- make_seurat_list(c("Sample_A", "Sample_B"))
    names(obj_list) <- c("Sample_A", "Sample_A")

    expect_error(mergeSeurat(obj_list), "duplicate names:\\s+Sample_A")
})

# ---- per-object preconditions --------------------------------------------

test_that("a zero-cell object is rejected", {
    obj_list <- list(Sample_A     = make_seurat("Sample_A"),
                     Sample_Empty = make_empty_seurat("Sample_Empty"))

    expect_error(mergeSeurat(obj_list), "zero cells:\\s+Sample_Empty")
})

test_that("mismatched default assays are rejected", {
    obj_list <- list(Sample_A = make_seurat("Sample_A"),
                     Sample_B = make_seurat("Sample_B", assay = "ADT"))

    expect_error(mergeSeurat(obj_list), "mismatched default assays")
})

test_that("a shared non-default assay name is accepted", {
    # the check is agreement, not RNA specifically
    obj_list <- list(Sample_A = make_seurat("Sample_A", assay = "ADT"),
                     Sample_B = make_seurat("Sample_B", assay = "ADT"))

    expect_output(res <- mergeSeurat(obj_list))
    expect_equal(Seurat::DefaultAssay(res), "ADT")
})

# ---- the merge itself ----------------------------------------------------

test_that("two objects merge into one, with barcodes namespaced by sample", {
    obj_list <- make_seurat_list(c("Sample_A", "Sample_B"), n_cells = 10)

    expect_output(res <- mergeSeurat(obj_list), "An object of class Seurat")

    expect_s4_class(res, "Seurat")
    expect_equal(ncol(res), 20)
    expect_setequal(colnames(res),
                    c(paste0("Sample_A_Cell", 1:10),
                      paste0("Sample_B_Cell", 1:10)))
})

test_that("Ident survives the merge and still identifies the sample", {
    obj_list <- make_seurat_list(c("Sample_A", "Sample_B"), n_cells = 10)

    expect_output(res <- mergeSeurat(obj_list))

    expect_setequal(unique(res$Ident), c("Sample_A", "Sample_B"))
    expect_equal(unname(table(res$Ident)[["Sample_A"]]), 10L)
})

test_that("colliding barcodes across samples are disambiguated, not dropped", {
    # every fixture object uses Cell1..Cell10 -- without add.cell.id these
    # would collide on merge
    obj_list <- make_seurat_list(c("Sample_A", "Sample_B", "Sample_C"),
                                 n_cells = 10)

    expect_output(res <- mergeSeurat(obj_list))

    expect_equal(ncol(res), 30)
    expect_equal(anyDuplicated(colnames(res)), 0L)
})

test_that("the merged object is returned invisibly", {
    obj_list <- make_seurat_list(c("Sample_A", "Sample_B"))

    capture.output(expect_invisible(mergeSeurat(obj_list)))
})

# ---- post-merge guards ---------------------------------------------------

test_that("objects with no Ident column are rejected before merging", {
    # Checked up front so the message is the package's own, and the same under
    # Seurat 4 and 5 -- left to the post-merge check, SeuratObject raised its
    # own error, worded differently by each version.
    obj_list <- list(Sample_A = make_seurat(NULL), Sample_B = make_seurat("Sample_B"))

    expect_error(mergeSeurat(obj_list), "no Ident column:\\s+Sample_A")
})

test_that("a sample whose Ident collapses into another is named", {
    a <- make_seurat("Sample_A")
    b <- make_seurat("Sample_B")
    b$Ident <- "Sample_A"          # b now claims to be A

    obj_list <- list(Sample_A = a, Sample_B = b)

    expect_error(mergeSeurat(obj_list),
                 "were not merged correctly:\\s+Sample_B")
})

test_that("a merge that loses cells is caught", {
    # merge() cannot realistically drop cells, so it is mocked -- the point is
    # that the guard fires rather than returning a short object.
    obj_list <- make_seurat_list(c("Sample_A", "Sample_B"), n_cells = 10)

    truncated <- make_seurat("Sample_A", n_cells = 5)
    truncated$Ident <- rep(c("Sample_A", "Sample_B"), length.out = 5)

    local_mocked_bindings(merge = function(...) truncated, .package = "base")

    expect_error(mergeSeurat(obj_list), "Cell count mismatch")
    expect_error(mergeSeurat(obj_list), "expected 20, got 5")
})

# ---- integration with seuratFrom10x() ------------------------------------

test_that("the list seuratFrom10x() returns merges without adjustment", {
    root <- local_root()
    local_mock_read10x()
    outs <- make_multi_run(root, "Pool1", c("Sample_A", "Sample_B"))

    seurat_list <- seuratFrom10x(path = outs, multi = TRUE,
                                 min.cells = 0, min.features = 0)

    expect_output(merged <- mergeSeurat(seurat_list))

    expect_equal(ncol(merged), 42)   # 21 cells per demuxed sample
    expect_setequal(unique(merged$Ident), c("Sample_A", "Sample_B"))
    expect_true(all(grepl("^Sample_[AB]_", colnames(merged))))
})
