# Internal QC helpers. Tests run inside the package namespace, so the @noRd
# functions are reachable directly.
#
# These three carry the argument surface every QC entry point shares, so they
# are tested once here rather than through each caller.
#
# Regexes use \\s+ where a message could wrap, since cli hard-wraps at the
# console width.

# ---- resolveOrganism -----------------------------------------------------

test_that("every supported organism resolves to itself", {
    for (code in c("hsa", "mmu", "rno", "dre")) {
        expect_equal(resolveOrganism(code), code)
    }
})

test_that("the full default vector resolves to the first element", {
    expect_equal(resolveOrganism(c("hsa", "mmu", "rno", "dre")), "hsa")
})

test_that("a unique prefix is accepted", {
    expect_equal(resolveOrganism("mm"), "mmu")
    expect_equal(resolveOrganism("d"), "dre")
})

test_that("an unsupported code is rejected and the supported set named", {
    expect_error(resolveOrganism("human"), "not supported")
    expect_error(resolveOrganism("human"), "hsa")
    expect_error(resolveOrganism("human"), "dre")
})

test_that("more than one organism is rejected", {
    expect_error(resolveOrganism(c("mmu", "hsa")),
                 "single organism code, not 2")
})

test_that("non-character, empty and NA input are rejected", {
    expect_error(resolveOrganism(1), "single organism code")
    expect_error(resolveOrganism(character(0)), "single organism code")
    expect_error(resolveOrganism(NA_character_), "single organism code")
})

test_that("errors are blamed on the caller, not on the helper", {
    caller <- function(organism) resolveOrganism(organism, call = environment())
    expect_error(caller("human"), "not supported")
})

# ---- genePatterns --------------------------------------------------------

test_that("each organism returns both patterns", {
    for (code in c("hsa", "mmu", "rno", "dre")) {
        p <- genePatterns(code)
        expect_named(p, c("mito", "ribo"))
        expect_type(p$mito, "character")
        expect_type(p$ribo, "character")
    }
})

test_that("patterns match the symbol convention they are for", {
    # the point of having four of these is that symbol case differs by species
    expect_true(grepl(genePatterns("hsa")$mito, "MT-CO1"))
    expect_false(grepl(genePatterns("hsa")$mito, "mt-Co1"))

    expect_true(grepl(genePatterns("mmu")$mito, "mt-Co1"))
    expect_true(grepl(genePatterns("mmu")$ribo, "Rps3"))
    expect_true(grepl(genePatterns("mmu")$ribo, "Rpl4"))

    expect_true(grepl(genePatterns("rno")$mito, "Mt-Co1"))
    expect_true(grepl(genePatterns("dre")$mito, "mt-co1"))
    expect_true(grepl(genePatterns("dre")$ribo, "rps3"))
})

test_that("the mouse mitochondrial pattern is hyphen-agnostic", {
    # GTF versions disagree on whether the symbol is mt-Co1 or mtCo1
    expect_true(grepl(genePatterns("mmu")$mito, "mt-Co1"))
    expect_true(grepl(genePatterns("mmu")$mito, "mtCo1"))
})

test_that("genePatterns defaults to human", {
    expect_identical(genePatterns(), genePatterns("hsa"))
})

test_that("genePatterns rejects what resolveOrganism rejects", {
    expect_error(genePatterns("rat"), "not supported")
})

# ---- validateQCObject ----------------------------------------------------

test_that("a valid object passes and returns invisibly", {
    obj <- make_seurat("Sample_A")
    expect_true(withVisible(validateQCObject(obj, "obj"))$value)
    expect_false(withVisible(validateQCObject(obj, "obj"))$visible)
})

test_that("a list of Seurat objects gets its own message", {
    obj_list <- make_seurat_list()
    expect_error(validateQCObject(obj_list, "obj_list"),
                 "is a list, not a single Seurat object")
})

test_that("the list message names the function the user actually called", {
    obj_list <- make_seurat_list()
    expect_error(validateQCObject(obj_list, "x", fn = "flagQC"),
                 "lapply\\(x, flagQC")
})

test_that("a data frame is not mistaken for a list of objects", {
    # data.frames are lists, so the list branch has to exclude them or the
    # hint sends the user somewhere useless
    expect_error(validateQCObject(data.frame(x = 1), "df"),
                 "must be a Seurat object")
})

test_that("non-Seurat input is rejected and its class named", {
    expect_error(validateQCObject(42, "x"), "not <numeric>")
    expect_error(validateQCObject("a", "x"), "not <character>")
})

test_that("an object that inherits from Seurat but fails validObject is rejected", {
    fake <- make_fake_seurat()
    expect_true(inherits(fake, "Seurat"))
    expect_error(validateQCObject(fake, "fake"), "not a structurally valid")
})

test_that("a zero-cell object is rejected", {
    expect_error(validateQCObject(make_empty_seurat(), "obj"), "zero cells")
})

test_that("errors name the argument as the caller wrote it", {
    expect_error(validateQCObject(42, "my_object"), "my_object")
})
