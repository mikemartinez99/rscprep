# runCellCycle() and the helpers behind it.
#
# The ortholog conversion is mocked throughout (see local_mock_orthologs()) so
# the tests are offline and deterministic -- the default `gprofiler` backend is
# a network call, and even the local backends are slow enough to notice. One
# test at the end runs the real conversion, skipped when orthogene is absent.
#
# Regexes use \\s+ where a message could wrap, since cli hard-wraps at the
# console width.

# ---- looksUnnormalised ---------------------------------------------------

test_that("raw counts are recognised as unnormalised", {
    expect_true(looksUnnormalised(make_cellcycle_seurat(normalise = FALSE)))
})

test_that("normalised data is not", {
    expect_false(looksUnnormalised(make_cellcycle_seurat(normalise = TRUE)))
})

test_that("an undecidable object returns NA rather than erroring", {
    # no non-zero values to inspect, so the heuristic has nothing to go on
    obj <- make_cellcycle_seurat(markers = "MCM5", n_bg = 5, n_cells = 4,
                                 normalise = FALSE)
    obj <- Seurat::SetAssayData(obj, slot = "data",
                                new.data = obj[["RNA"]]@data * 0)

    expect_true(is.na(looksUnnormalised(obj)))
})

# ---- cellCycleFeatures ---------------------------------------------------

test_that("human needs no conversion at all", {
    skip_if_not_installed("orthogene")
    # mock set to explode: reaching it would be the bug
    testthat::local_mocked_bindings(
        convert_orthologs = function(...) stop("must not be called for human"),
        .package = "orthogene"
    )

    res <- cellCycleFeatures("hsa", "gprofiler")

    expect_identical(res$s, Seurat::cc.genes$s.genes)
    expect_identical(res$g2m, Seurat::cc.genes$g2m.genes)
})

test_that("a non-human organism goes through conversion", {
    skip_if_not_installed("orthogene")
    local_mock_orthologs()

    res <- cellCycleFeatures("mmu", "gprofiler")

    expect_equal(res$s, pseudo_ortholog(Seurat::cc.genes$s.genes))
    expect_equal(res$g2m, pseudo_ortholog(Seurat::cc.genes$g2m.genes))
})

# ---- cellCycleOrthologs --------------------------------------------------

test_that("both marker sets are converted in a single call", {
    # each call is a network round trip on the default backend, and the sets
    # do not overlap, so one request should do the work of two
    skip_if_not_installed("orthogene")
    calls <- local_mock_orthologs()

    cellCycleOrthologs(c("AAA", "BBB"), c("CCC"), "mmu", "gprofiler")

    expect_length(calls$gene_df, 1)
    expect_setequal(calls$gene_df[[1]], c("AAA", "BBB", "CCC"))
})

test_that("results are matched by input symbol, not by position", {
    # the mock returns its rows reversed; relying on order would swap the sets
    skip_if_not_installed("orthogene")
    local_mock_orthologs(shuffle = TRUE)

    res <- cellCycleOrthologs(c("AAA", "BBB"), c("CCC", "DDD"), "mmu", "gprofiler")

    expect_equal(res$s, pseudo_ortholog(c("AAA", "BBB")))
    expect_equal(res$g2m, pseudo_ortholog(c("CCC", "DDD")))
})

test_that("genes with no one-to-one ortholog are dropped", {
    skip_if_not_installed("orthogene")
    local_mock_orthologs(drop = c("BBB", "DDD"))

    res <- cellCycleOrthologs(c("AAA", "BBB"), c("CCC", "DDD"), "mmu", "gprofiler")

    expect_equal(res$s, pseudo_ortholog("AAA"))
    expect_equal(res$g2m, pseudo_ortholog("CCC"))
})

test_that("dropped genes are reported", {
    skip_if_not_installed("orthogene")
    local_mock_orthologs(drop = "BBB")

    expect_message(
        cellCycleOrthologs(c("AAA", "BBB"), "CCC", "mmu", "gprofiler"),
        "1 marker had no one-to-one"
    )
})

test_that("nothing is reported when every gene converts", {
    skip_if_not_installed("orthogene")
    local_mock_orthologs()

    expect_no_message(
        cellCycleOrthologs("AAA", "BBB", "mmu", "gprofiler"),
        message = "had no one-to-one"
    )
})

test_that("an empty S set after conversion is an error", {
    skip_if_not_installed("orthogene")
    local_mock_orthologs(drop = c("AAA", "BBB"))

    expect_error(cellCycleOrthologs(c("AAA", "BBB"), "CCC", "mmu", "gprofiler"),
                 "No \"mouse\" orthologs were found for the\\s+S markers")
})

test_that("an empty G2/M set after conversion is an error", {
    skip_if_not_installed("orthogene")
    local_mock_orthologs(drop = "CCC")

    expect_error(cellCycleOrthologs("AAA", "CCC", "mmu", "gprofiler"),
                 "G2/M markers")
})

test_that("a conversion failure is wrapped with the offline hint", {
    skip_if_not_installed("orthogene")
    testthat::local_mocked_bindings(
        convert_orthologs = function(...) stop("could not reach gprofiler"),
        .package = "orthogene"
    )

    expect_error(cellCycleOrthologs("AAA", "BBB", "mmu", "gprofiler"),
                 "Ortholog conversion to \"mouse\" failed")
    expect_error(cellCycleOrthologs("AAA", "BBB", "mmu", "gprofiler"),
                 "homologene")
})

test_that("each organism is sent to orthogene under the species name it expects", {
    skip_if_not_installed("orthogene")
    calls <- local_mock_orthologs()

    for (code in c("mmu", "rno", "dre")) {
        cellCycleOrthologs("AAA", "BBB", code, "gprofiler")
    }

    expect_equal(calls$output_species, c("mouse", "rat", "zebrafish"))
})

# ---- runCellCycle: what it adds ------------------------------------------

test_that("all three columns are added", {
    res <- runCellCycle(make_cellcycle_seurat(), "hsa")

    expect_true(all(c("S.Score", "G2M.Score", "Phase") %in%
                        colnames(res@meta.data)))
})

test_that("Phase is one of the three calls, for every cell", {
    res <- runCellCycle(make_cellcycle_seurat(), "hsa")

    expect_true(all(res$Phase %in% c("G1", "S", "G2M")))
    expect_equal(sum(table(res$Phase)), ncol(res))
})

test_that("nothing is filtered and identities survive", {
    obj <- make_cellcycle_seurat()
    res <- runCellCycle(obj, "hsa")

    expect_s4_class(res, "Seurat")
    expect_equal(ncol(res), ncol(obj))
    expect_equal(colnames(res), colnames(obj))
    expect_equal(Seurat::Idents(res), Seurat::Idents(obj))
})

test_that("the scores are numeric", {
    res <- runCellCycle(make_cellcycle_seurat(), "hsa")

    expect_type(res$S.Score, "double")
    expect_type(res$G2M.Score, "double")
})

# ---- runCellCycle: normalisation -----------------------------------------

test_that("a raw object is warned about and then normalised", {
    obj <- make_cellcycle_seurat(normalise = FALSE)

    expect_warning(runCellCycle(obj, "hsa"), "still holds raw counts")
    res <- suppressWarnings(runCellCycle(obj, "hsa"))

    expect_false(looksUnnormalised(res))
})

test_that("the normalisation the object gets is the ordinary one", {
    obj <- make_cellcycle_seurat(normalise = FALSE)

    res <- suppressWarnings(runCellCycle(obj, "hsa"))

    expect_equal(Seurat::GetAssayData(res, slot = "data"),
                 Seurat::GetAssayData(Seurat::NormalizeData(obj, verbose = FALSE),
                                      slot = "data"))
})

test_that("normalising leaves the counts layer alone", {
    obj <- make_cellcycle_seurat(normalise = FALSE)

    res <- suppressWarnings(runCellCycle(obj, "hsa"))

    expect_equal(Seurat::GetAssayData(res, slot = "counts"),
                 Seurat::GetAssayData(obj, slot = "counts"))
})

test_that("an already-normalised object is left exactly as it was", {
    obj <- make_cellcycle_seurat(normalise = TRUE)

    res <- expect_no_warning(runCellCycle(obj, "hsa"))

    expect_equal(Seurat::GetAssayData(res, slot = "data"),
                 Seurat::GetAssayData(obj, slot = "data"))
})

# ---- runCellCycle: markers against the object ----------------------------

test_that("an object with no markers at all is an error", {
    obj <- make_cellcycle_seurat(markers = character(0))

    expect_error(runCellCycle(obj, "hsa"), "None of the \"hsa\"")
    expect_error(runCellCycle(obj, "hsa"), "Matched 0 of 43 S")
})

test_that("the error counts against the full marker set, not the matched one", {
    obj <- make_cellcycle_seurat(markers = character(0))

    expect_error(runCellCycle(obj, "hsa"), "0 of 54 G2/M")
})

test_that("few matched markers is a warning, not an error", {
    obj <- make_cellcycle_seurat(
        markers = c(Seurat::cc.genes$s.genes[1:4], Seurat::cc.genes$g2m.genes[1:4])
    )

    expect_warning(runCellCycle(obj, "hsa"), "Few markers matched: 4 S, 4 G2/M")
    expect_true("Phase" %in% colnames(suppressWarnings(runCellCycle(obj, "hsa"))@meta.data))
})

test_that("a comfortable overlap warns about nothing", {
    expect_no_warning(runCellCycle(make_cellcycle_seurat(), "hsa"))
})

test_that("a scoring failure is wrapped with a cause worth reading", {
    # too small a background pool for the control set -- Seurat's own message
    # is about sampling and says nothing about cell cycle
    obj <- make_cellcycle_seurat(n_bg = 40, n_cells = 20)

    expect_error(runCellCycle(obj, "hsa"), "CellCycleScoring.*failed")
    expect_error(runCellCycle(obj, "hsa"), "too few genes or cells")
})

# ---- runCellCycle: validation --------------------------------------------

test_that("a list of objects is rejected with a usable hint", {
    obj_list <- make_seurat_list()

    expect_error(runCellCycle(obj_list), "is a list, not a single Seurat object")
    expect_error(runCellCycle(obj_list), "lapply\\(obj_list, runCellCycle")
})

test_that("non-Seurat input is rejected", {
    expect_error(runCellCycle(42), "must be a Seurat object")
})

test_that("errors name the argument as the caller wrote it", {
    my_object <- 42
    expect_error(runCellCycle(my_object), "my_object")
})

test_that("an unsupported organism is rejected", {
    expect_error(runCellCycle(make_cellcycle_seurat(), "human"), "not supported")
})

test_that("an unsupported method is rejected", {
    expect_error(runCellCycle(make_cellcycle_seurat(), "hsa", method = "nope"),
                 "'arg' should be one of")
})

test_that("organism defaults to human", {
    obj <- make_cellcycle_seurat()

    expect_equal(runCellCycle(obj)$Phase, runCellCycle(obj, "hsa")$Phase)
})

# ---- runCellCycle: the non-human path ------------------------------------

test_that("a non-human organism scores against converted markers", {
    skip_if_not_installed("orthogene")
    local_mock_orthologs()

    obj <- make_cellcycle_seurat(
        markers = pseudo_ortholog(c(Seurat::cc.genes$s.genes,
                                    Seurat::cc.genes$g2m.genes))
    )

    res <- runCellCycle(obj, "mmu")

    expect_true(all(res$Phase %in% c("G1", "S", "G2M")))
})

test_that("human markers in a mouse-coded call find nothing", {
    skip_if_not_installed("orthogene")
    local_mock_orthologs()

    # object carries human symbols, but conversion produced mouse ones
    obj <- make_cellcycle_seurat()

    expect_error(runCellCycle(obj, "mmu"), "None of the \"mmu\"")
})

# ---- the real conversion -------------------------------------------------

test_that("the real ortholog conversion works end to end", {
    skip_if_not_installed("orthogene")
    skip_on_cran()

    # homologene resolves from a bundled mapping, so this needs no network
    res <- cellCycleOrthologs(Seurat::cc.genes$s.genes,
                              Seurat::cc.genes$g2m.genes,
                              "mmu", "homologene")

    expect_type(res$s, "character")
    expect_type(res$g2m, "character")
    # mouse symbols are title case where human ones are upper
    expect_true(all(res$s != toupper(res$s)))
    # conversion drops some genes but should not gut the set
    expect_gt(length(res$s), 30)
    expect_gt(length(res$g2m), 40)
    expect_length(intersect(res$s, res$g2m), 0)
})
