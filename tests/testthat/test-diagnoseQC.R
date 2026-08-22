# diagnoseQC() adds the three per-cell QC metrics. Fixtures use all-1 counts,
# which makes every expected value arithmetic rather than approximate: with
# one count per gene per cell, pctMito is 100 * n_mito / n_genes exactly, and
# log10GenesPerUMI is log10(n)/log10(n) = 1.
#
# Regexes use \\s+ where a message could wrap, since cli hard-wraps at the
# console width.

# ---- the metrics ---------------------------------------------------------

test_that("all three metrics are added", {
    obj <- make_pattern_seurat()

    res <- expect_no_warning(diagnoseQC(obj, "mmu"))

    expect_true(all(c("pctMito", "pctRibo", "log10GenesPerUMI") %in%
                        colnames(res@meta.data)))
})

test_that("percentages are the exact share of counts the patterns match", {
    # 3 mito + 2 ribo + 15 other = 20 genes, one count each
    obj <- make_pattern_seurat(mito = c("mt-Co1", "mt-Co2", "mt-Nd1"),
                               ribo = c("Rps3", "Rpl4"),
                               n_other = 15)

    res <- diagnoseQC(obj, "mmu")

    expect_true(all(res$pctMito == 15))
    expect_true(all(res$pctRibo == 10))
})

test_that("complexity is log10(nFeature) / log10(nCount)", {
    obj <- make_pattern_seurat()
    res <- diagnoseQC(obj, "mmu")

    expect_equal(res$log10GenesPerUMI,
                 log10(res$nFeature_RNA) / log10(res$nCount_RNA))
    # every gene detected exactly once, so the ratio is 1
    expect_true(all(res$log10GenesPerUMI == 1))
})

test_that("nothing is filtered", {
    obj <- make_pattern_seurat(n_cells = 20)
    res <- diagnoseQC(obj, "mmu")

    expect_equal(ncol(res), 20)
    expect_equal(nrow(res), nrow(obj))
    expect_equal(colnames(res), colnames(obj))
})

test_that("a Seurat object is returned", {
    expect_s4_class(diagnoseQC(make_pattern_seurat(), "mmu"), "Seurat")
})

test_that("re-running overwrites rather than duplicating", {
    obj <- diagnoseQC(make_pattern_seurat(), "mmu")
    res <- diagnoseQC(obj, "mmu")

    expect_equal(sum(colnames(res@meta.data) == "pctMito"), 1)
    expect_equal(res$pctMito, obj$pctMito)
})

test_that("re-running with a different organism replaces the values", {
    obj <- diagnoseQC(make_pattern_seurat(), "mmu")
    expect_true(all(obj$pctMito == 15))

    res <- suppressWarnings(diagnoseQC(obj, "hsa"))

    # human patterns match no mouse symbol, so the columns reset to zero
    expect_true(all(res$pctMito == 0))
})

# ---- organism resolution -------------------------------------------------

test_that("organism defaults to human", {
    obj <- make_pattern_seurat(mito = "MT-CO1", ribo = "RPS3", n_other = 18)
    res <- diagnoseQC(obj)

    expect_true(all(res$pctMito == 5))
})

test_that("a unique organism prefix is accepted", {
    obj <- make_pattern_seurat()
    expect_equal(diagnoseQC(obj, "mm")$pctMito, diagnoseQC(obj, "mmu")$pctMito)
})

test_that("an unsupported organism is rejected before any work is done", {
    expect_error(diagnoseQC(make_pattern_seurat(), "human"), "not supported")
})

test_that("the reported organism is the resolved one, not the default vector", {
    # organism defaults to a 4-element vector; reporting it unresolved would
    # name all four as the one in use
    obj <- make_pattern_seurat(mito = character(0), ribo = character(0))

    w <- testthat::capture_warnings(diagnoseQC(obj))
    expect_true(all(grepl("organism \"hsa\"", w)))
    expect_false(any(grepl("mmu|rno|dre", w)))
})

# ---- patterns that match nothing -----------------------------------------

test_that("a pattern matching no feature warns and names the pattern", {
    obj <- make_pattern_seurat()  # mouse symbols

    w <- testthat::capture_warnings(diagnoseQC(obj, "hsa"))
    expect_length(w, 2)
    expect_true(all(grepl("matched 0 features", w)))
    expect_true(any(grepl("\\^MT-", w)))
    expect_true(any(grepl("\\^RP\\[SL\\]", w)))
})

test_that("the warning names the column that will be empty", {
    obj <- make_pattern_seurat(mito = character(0), n_other = 18)

    expect_warning(diagnoseQC(obj, "mmu"), "pctMito will be 0")
})

test_that("a zero-match pattern still produces a column, of zeroes", {
    obj <- make_pattern_seurat(mito = character(0), n_other = 18)

    res <- suppressWarnings(diagnoseQC(obj, "mmu"))

    expect_true("pctMito" %in% colnames(res@meta.data))
    expect_true(all(res$pctMito == 0))
})

test_that("one failing pattern does not stop the other from being calculated", {
    obj <- make_pattern_seurat(mito = character(0), ribo = c("Rps3", "Rpl4"),
                               n_other = 18)

    res <- suppressWarnings(diagnoseQC(obj, "mmu"))

    expect_true(all(res$pctMito == 0))
    expect_true(all(res$pctRibo == 10))
})

test_that("a matching pattern warns about neither set", {
    expect_no_warning(diagnoseQC(make_pattern_seurat(), "mmu"))
})

# ---- the count columns ---------------------------------------------------

test_that("count columns are read from the default assay", {
    obj <- make_pattern_seurat(assay = "SCT")

    expect_equal(Seurat::DefaultAssay(obj), "SCT")
    expect_true("nCount_SCT" %in% colnames(obj@meta.data))
    expect_false("nCount_RNA" %in% colnames(obj@meta.data))
    res <- diagnoseQC(obj, "mmu")

    expect_true("log10GenesPerUMI" %in% colnames(res@meta.data))
    expect_true(all(res$log10GenesPerUMI == 1))
})

test_that("missing count columns are a clear error, not a cryptic one", {
    obj <- make_pattern_seurat()
    obj@meta.data$nCount_RNA   <- NULL
    obj@meta.data$nFeature_RNA <- NULL

    expect_error(diagnoseQC(obj, "mmu"), "missing the count columns")
    expect_error(diagnoseQC(obj, "mmu"), "nCount_RNA")
})

test_that("the error names only the columns actually absent", {
    obj <- make_pattern_seurat()
    obj@meta.data$nCount_RNA <- NULL

    expect_error(diagnoseQC(obj, "mmu"), "count column \"nCount_RNA\"")
})

# ---- non-finite complexity -----------------------------------------------

test_that("a single-UMI cell warns instead of silently producing NaN", {
    # log10(1) is 0, so the complexity denominator vanishes
    genes <- c("mt-Co1", paste0("Gene", 1:9))
    m <- matrix(1, nrow = 10, ncol = 5,
                dimnames = list(genes, paste0("Cell", 1:5)))
    m[, 1]  <- 0
    m[1, 1] <- 1
    obj <- Seurat::CreateSeuratObject(Matrix::Matrix(m, sparse = TRUE),
                                      min.cells = 0, min.features = 0)

    w <- testthat::capture_warnings(diagnoseQC(obj, "mmu"))
    expect_true(any(grepl("non-finite", w)))
    res <- suppressWarnings(diagnoseQC(obj, "mmu"))

    expect_false(is.finite(res$log10GenesPerUMI[1]))
    expect_true(all(is.finite(res$log10GenesPerUMI[-1])))
})

# ---- validation ----------------------------------------------------------

test_that("a list of objects is rejected with a usable hint", {
    obj_list <- make_seurat_list()

    expect_error(diagnoseQC(obj_list), "is a list, not a single Seurat object")
    expect_error(diagnoseQC(obj_list), "lapply\\(obj_list, diagnoseQC")
})

test_that("non-Seurat input is rejected", {
    expect_error(diagnoseQC(42), "must be a Seurat object")
})

test_that("errors name the argument as the caller wrote it", {
    my_object <- 42
    expect_error(diagnoseQC(my_object), "my_object")
})
