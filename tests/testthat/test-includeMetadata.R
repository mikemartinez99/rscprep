# includeMetadata() joins a sample-level table onto the named list
# seuratFrom10x() returns. Its defining contract is completeness: every object
# must get every column. Partial annotation is rejected precisely because
# Seurat's merge() does NOT reject it -- it fills the gap with NA and says
# nothing (see the final section, which pins that behaviour).
#
# Regexes use \\s+ where a message could wrap, since cli hard-wraps at the
# console width.

# ---- sample key resolution -----------------------------------------------

test_that("samples can be keyed by a Sample column", {
    obj_list <- make_seurat_list(c("Sample_A", "Sample_B"))
    sheet    <- make_sample_sheet(c("Sample_A", "Sample_B"), key = "column")

    res <- includeMetadata(obj_list, sheet)

    expect_true(all(res$Sample_A$Treatment == "Control"))
    expect_true(all(res$Sample_B$Treatment == "Bleomycin"))
})

test_that("samples can be keyed by rownames", {
    obj_list <- make_seurat_list(c("Sample_A", "Sample_B"))
    sheet    <- make_sample_sheet(c("Sample_A", "Sample_B"), key = "rownames")

    res <- includeMetadata(obj_list, sheet)

    expect_true(all(res$Sample_A$Treatment == "Control"))
    expect_true(all(res$Sample_B$Treatment == "Bleomycin"))
})

test_that("a Sample column wins over rownames", {
    obj_list <- make_seurat_list(c("Sample_A", "Sample_B"))
    sheet    <- make_sample_sheet(c("Sample_A", "Sample_B"), key = "column")
    # rownames deliberately disagree with the column, and are ignored
    rownames(sheet) <- c("Sample_B", "Sample_A")

    res <- includeMetadata(obj_list, sheet)

    expect_true(all(res$Sample_A$Treatment == "Control"))
})

test_that("the Sample key column is not attached as metadata", {
    obj_list <- make_seurat_list(c("Sample_A", "Sample_B"))
    sheet    <- make_sample_sheet(c("Sample_A", "Sample_B"), key = "column")

    res <- includeMetadata(obj_list, sheet)

    expect_false("Sample" %in% colnames(res$Sample_A[[]]))
    expect_true("Treatment" %in% colnames(res$Sample_A[[]]))
})

test_that("automatic rownames are rejected rather than matched as \"1\", \"2\"", {
    obj_list <- make_seurat_list(c("Sample_A", "Sample_B"))
    sheet    <- data.frame(Treatment = c("Control", "Bleomycin"))

    expect_error(includeMetadata(obj_list, sheet),
                 "Can't tell which sample each row")
})

test_that("a tibble without a Sample column is rejected", {
    skip_if_not_installed("tibble")

    obj_list <- make_seurat_list(c("Sample_A", "Sample_B"))
    sheet    <- tibble::tibble(Treatment = c("Control", "Bleomycin"))

    expect_error(includeMetadata(obj_list, sheet),
                 "Can't tell which sample each row")
})

test_that("a tibble with a Sample column works", {
    skip_if_not_installed("tibble")

    obj_list <- make_seurat_list(c("Sample_A", "Sample_B"))
    sheet    <- tibble::tibble(Sample    = c("Sample_A", "Sample_B"),
                               Treatment = c("Control", "Bleomycin"))

    res <- includeMetadata(obj_list, sheet)

    expect_true(all(res$Sample_B$Treatment == "Bleomycin"))
})

# ---- the happy path ------------------------------------------------------

test_that("a value is broadcast across every cell of its own object", {
    obj_list <- make_seurat_list(c("Sample_A", "Sample_B"), n_cells = 10)
    sheet    <- make_sample_sheet(c("Sample_A", "Sample_B"))

    res <- includeMetadata(obj_list, sheet)

    expect_length(res$Sample_A$Treatment, 10)
    expect_equal(unique(res$Sample_A$Treatment), "Control")
    expect_equal(unique(res$Sample_B$Treatment), "Bleomycin")
})

test_that("several columns are attached at once", {
    obj_list <- make_seurat_list(c("Sample_A", "Sample_B"))
    sheet    <- make_sample_sheet(c("Sample_A", "Sample_B"),
                                  Treatment = c("Control", "Bleomycin"),
                                  Timepoint = c("D0", "D21"),
                                  Batch     = c("B1", "B1"))

    res <- includeMetadata(obj_list, sheet)

    expect_true(all(c("Treatment", "Timepoint", "Batch") %in%
                        colnames(res$Sample_A[[]])))
    expect_equal(unique(res$Sample_B$Timepoint), "D21")
})

test_that("list names and order are preserved", {
    obj_list <- make_seurat_list(c("Sample_C", "Sample_A", "Sample_B"))
    sheet    <- make_sample_sheet(c("Sample_A", "Sample_B", "Sample_C"),
                                  Treatment = c("x", "y", "z"))

    res <- includeMetadata(obj_list, sheet)

    expect_equal(names(res), c("Sample_C", "Sample_A", "Sample_B"))
    # and each object still got its OWN row, not the sheet's row order
    expect_equal(unique(res$Sample_C$Treatment), "z")
    expect_equal(unique(res$Sample_A$Treatment), "x")
})

test_that("column types are preserved", {
    obj_list <- make_seurat_list(c("Sample_A", "Sample_B"))
    sheet    <- make_sample_sheet(c("Sample_A", "Sample_B"),
                                  Treatment = c("Control", "Bleomycin"),
                                  Dose      = c(0, 2.5),
                                  Responder = c(TRUE, FALSE))

    res <- includeMetadata(obj_list, sheet)

    expect_type(res$Sample_A$Treatment, "character")
    expect_type(res$Sample_A$Dose, "double")
    expect_type(res$Sample_A$Responder, "logical")
    expect_equal(unique(res$Sample_B$Dose), 2.5)
})

test_that("factor levels are preserved on every object", {
    # each object only carries one of the levels, but keeps the full set, so
    # per-sample subsets stay comparable
    obj_list <- make_seurat_list(c("Sample_A", "Sample_B"))
    sheet    <- make_sample_sheet(
        c("Sample_A", "Sample_B"),
        Treatment = factor(c("Control", "Bleomycin"),
                           levels = c("Control", "Bleomycin"))
    )

    res <- includeMetadata(obj_list, sheet)

    expect_s3_class(res$Sample_A$Treatment, "factor")
    expect_equal(levels(res$Sample_A$Treatment), c("Control", "Bleomycin"))
    expect_equal(levels(res$Sample_B$Treatment), c("Control", "Bleomycin"))
})

test_that("the annotated list is returned, not modified in place", {
    obj_list <- make_seurat_list(c("Sample_A", "Sample_B"))
    sheet    <- make_sample_sheet(c("Sample_A", "Sample_B"))

    res <- includeMetadata(obj_list, sheet)

    expect_false("Treatment" %in% colnames(obj_list$Sample_A[[]]))
    expect_true("Treatment" %in% colnames(res$Sample_A[[]]))
})

# ---- completeness, part 1: every sample needs a row -----------------------

test_that("a sample with no row is an error naming it", {
    obj_list <- make_seurat_list(c("Sample_A", "Sample_B"))
    sheet    <- make_sample_sheet("Sample_A")

    expect_error(includeMetadata(obj_list, sheet), "needs a row")
    expect_error(includeMetadata(obj_list, sheet), "No row for:\\s+.Sample_B.")
})

test_that("the missing-row error fires even when other samples are complete", {
    obj_list <- make_seurat_list(c("Sample_A", "Sample_B", "Sample_C"))
    sheet    <- make_sample_sheet(c("Sample_A", "Sample_C"),
                                  Treatment = c("Control", "Bleomycin"))

    expect_error(includeMetadata(obj_list, sheet), "Sample_B")
})

# ---- completeness, part 2: no gaps within a column -----------------------

test_that("NA in a column errors, naming the column and the sample", {
    obj_list <- make_seurat_list(c("Sample_A", "Sample_B"))
    sheet    <- make_sample_sheet(c("Sample_A", "Sample_B"),
                                  Treatment = c("Control", NA))

    expect_error(includeMetadata(obj_list, sheet), "must be populated for every sample")
    expect_error(includeMetadata(obj_list, sheet), "Treatment.*missing for")
    expect_error(includeMetadata(obj_list, sheet), "Sample_B")
})

test_that("an empty string counts as missing", {
    # a blank CSV cell reads as "", not NA, and is just as unusable
    obj_list <- make_seurat_list(c("Sample_A", "Sample_B"))
    sheet    <- make_sample_sheet(c("Sample_A", "Sample_B"),
                                  Treatment = c("Control", ""))

    expect_error(includeMetadata(obj_list, sheet), "must be populated for every sample")
})

test_that("a whitespace-only string counts as missing", {
    obj_list <- make_seurat_list(c("Sample_A", "Sample_B"))
    sheet    <- make_sample_sheet(c("Sample_A", "Sample_B"),
                                  Treatment = c("Control", "   "))

    expect_error(includeMetadata(obj_list, sheet), "must be populated for every sample")
})

test_that("NA in a numeric or logical column is caught too", {
    obj_list <- make_seurat_list(c("Sample_A", "Sample_B"))

    numeric_sheet <- make_sample_sheet(c("Sample_A", "Sample_B"),
                                       Dose = c(1.5, NA))
    logical_sheet <- make_sample_sheet(c("Sample_A", "Sample_B"),
                                       Responder = c(TRUE, NA))

    expect_error(includeMetadata(obj_list, numeric_sheet), "Dose")
    expect_error(includeMetadata(obj_list, logical_sheet), "Responder")
})

test_that("every incomplete column is reported in one message", {
    # fixing a 40-sample sheet one error per run is miserable
    obj_list <- make_seurat_list(c("Sample_A", "Sample_B"))
    sheet    <- make_sample_sheet(c("Sample_A", "Sample_B"),
                                  Treatment = c("Control", NA),
                                  Timepoint = c(NA, "D21"),
                                  Batch     = c("B1", "B1"))

    err <- expect_error(includeMetadata(obj_list, sheet))
    msg <- conditionMessage(err)

    expect_match(msg, "Treatment")
    expect_match(msg, "Timepoint")
    expect_no_match(msg, "Batch")
})

test_that("blanks are only checked for samples that are actually loaded", {
    # an incomplete row for a sample we did not load is irrelevant
    obj_list <- make_seurat_list("Sample_A")
    sheet    <- make_sample_sheet(c("Sample_A", "Sample_B"),
                                  Treatment = c("Control", NA))

    expect_warning(res <- includeMetadata(obj_list, sheet), "Ignoring")
    expect_equal(unique(res$Sample_A$Treatment), "Control")
})

# ---- extra rows ----------------------------------------------------------

test_that("rows matching no object warn and are ignored", {
    # the master-sheet workflow: one sheet for the experiment, a subset loaded
    obj_list <- make_seurat_list(c("Sample_A", "Sample_B"))
    sheet    <- make_sample_sheet(c("Sample_A", "Sample_B", "Sample_C", "Sample_D"),
                                  Treatment = c("w", "x", "y", "z"))

    expect_warning(res <- includeMetadata(obj_list, sheet),
                   "Ignoring 2 rows")
    expect_length(res, 2)
    expect_equal(unique(res$Sample_A$Treatment), "w")
    expect_equal(unique(res$Sample_B$Treatment), "x")
})

test_that("a mistyped sample name is caught despite extra rows being tolerated", {
    obj_list <- make_seurat_list(c("Sample_A", "Sample_B"))
    sheet    <- make_sample_sheet(c("Sample_A", "Sanple_B"),
                                  Treatment = c("Control", "Bleomycin"))

    # the typo shows up as a sample with no row -- a hard error
    expect_error(suppressWarnings(includeMetadata(obj_list, sheet)), "Sample_B")
})

# ---- collisions ----------------------------------------------------------

test_that("an existing column is refused by default", {
    obj_list <- make_seurat_list(c("Sample_A", "Sample_B"))
    sheet    <- make_sample_sheet(c("Sample_A", "Sample_B"))
    obj_list <- includeMetadata(obj_list, sheet)

    expect_error(includeMetadata(obj_list, sheet), "already present")
    expect_error(includeMetadata(obj_list, sheet), "overwrite = TRUE")
})

test_that("overwrite = TRUE replaces an existing column", {
    obj_list <- make_seurat_list(c("Sample_A", "Sample_B"))
    obj_list <- includeMetadata(obj_list, make_sample_sheet(c("Sample_A", "Sample_B")))

    corrected <- make_sample_sheet(c("Sample_A", "Sample_B"),
                                   Treatment = c("Vehicle", "Bleomycin"))
    res <- includeMetadata(obj_list, corrected, overwrite = TRUE)

    expect_equal(unique(res$Sample_A$Treatment), "Vehicle")
})

test_that("a collision with a Seurat-supplied column is caught", {
    obj_list <- make_seurat_list(c("Sample_A", "Sample_B"))
    sheet    <- make_sample_sheet(c("Sample_A", "Sample_B"),
                                  orig.ident = c("x", "y"))

    expect_error(includeMetadata(obj_list, sheet), "already present")
})

test_that("Ident is reserved even under overwrite = TRUE", {
    obj_list <- make_seurat_list(c("Sample_A", "Sample_B"))
    sheet    <- make_sample_sheet(c("Sample_A", "Sample_B"),
                                  Ident = c("x", "y"))

    expect_error(includeMetadata(obj_list, sheet, overwrite = TRUE), "reserved")
})

# ---- input validation ----------------------------------------------------

test_that("objList must be a non-empty named list of Seurat objects", {
    sheet <- make_sample_sheet(c("Sample_A", "Sample_B"))

    expect_error(includeMetadata(list(), sheet), "is empty")
    expect_error(includeMetadata(unname(make_seurat_list()), sheet),
                 "must be fully named")
    expect_error(
        includeMetadata(list(Sample_A = make_seurat("Sample_A"),
                             Sample_B = "not a seurat object"), sheet),
        "not a .*Seurat.* object"
    )
})

test_that("a single Seurat object is rejected with a pointed message", {
    expect_error(
        includeMetadata(make_seurat("Sample_A"),
                        make_sample_sheet("Sample_A")),
        "not a single object"
    )
})

test_that("metadata must be a data frame with rows and columns", {
    obj_list <- make_seurat_list(c("Sample_A", "Sample_B"))

    expect_error(includeMetadata(obj_list, "not a data frame"),
                 "must be a .*data.frame")
    expect_error(
        includeMetadata(obj_list, data.frame(Sample = character(0),
                                             Treatment = character(0))),
        "has no rows"
    )
    expect_error(
        includeMetadata(obj_list, data.frame(Sample = c("Sample_A", "Sample_B"))),
        "no metadata columns"
    )
})

test_that("duplicate sample names in the sheet are rejected", {
    obj_list <- make_seurat_list(c("Sample_A", "Sample_B"))
    sheet    <- make_sample_sheet(c("Sample_A", "Sample_A"),
                                  Treatment = c("Control", "Bleomycin"))

    expect_error(includeMetadata(obj_list, sheet),
                 "duplicate sample names:\\s+.Sample_A.")
})

test_that("blank sample names in the sheet are rejected", {
    obj_list <- make_seurat_list(c("Sample_A", "Sample_B"))
    sheet    <- make_sample_sheet(c("Sample_A", ""),
                                  Treatment = c("Control", "Bleomycin"))

    expect_error(includeMetadata(obj_list, sheet), "missing or blank sample names")
})

# ---- column names --------------------------------------------------------

test_that("syntactically invalid column names are refused, not silently renamed", {
    # Seurat runs make.names() on metadata column names without asking, so
    # "Cell Type" would land as `Cell.Type` and stop matching the sheet header
    obj_list <- make_seurat_list(c("Sample_A", "Sample_B"))
    sheet    <- data.frame(Sample = c("Sample_A", "Sample_B"),
                           check.names = FALSE)
    sheet[["Cell Type"]] <- c("T cell", "B cell")

    expect_error(includeMetadata(obj_list, sheet), "not syntactically valid")
    expect_error(includeMetadata(obj_list, sheet), "Cell.Type")
})

test_that("a column name starting with a digit is refused", {
    obj_list <- make_seurat_list(c("Sample_A", "Sample_B"))
    sheet    <- data.frame(Sample = c("Sample_A", "Sample_B"),
                           check.names = FALSE)
    sheet[["2Group"]] <- c("a", "b")

    expect_error(includeMetadata(obj_list, sheet), "not syntactically valid")
})

test_that("valid column names with dots and underscores are accepted", {
    obj_list <- make_seurat_list(c("Sample_A", "Sample_B"))
    sheet    <- make_sample_sheet(c("Sample_A", "Sample_B"),
                                  cell.type   = c("T", "B"),
                                  donor_id    = c("D1", "D2"))

    res <- includeMetadata(obj_list, sheet)

    expect_true(all(c("cell.type", "donor_id") %in% colnames(res$Sample_A[[]])))
})

test_that("braces in names do not blow up the error path", {
    # cli_abort() glue-interpolates the bullets it is handed, so an unescaped
    # brace in a sample name would be evaluated as an expression and replace
    # the real message with a cli internal error
    obj_list <- make_seurat_list(c("S{A}", "Sample_B"))
    sheet    <- data.frame(Sample    = c("S{A}", "Sample_B"),
                           Treatment = c(NA, "Bleomycin"))

    expect_error(includeMetadata(obj_list, sheet),
                 "must be populated for every sample")
    expect_error(includeMetadata(obj_list, sheet), "S\\{A\\}")
})

test_that("braces in a metadata value are attached verbatim", {
    obj_list <- make_seurat_list(c("Sample_A", "Sample_B"))
    sheet    <- make_sample_sheet(c("Sample_A", "Sample_B"),
                                  Treatment = c("Drug{1}", "Vehicle"))

    res <- includeMetadata(obj_list, sheet)

    expect_equal(unique(res$Sample_A$Treatment), "Drug{1}")
})

# ---- smaller edges -------------------------------------------------------

test_that("a one-sample list and a one-row sheet work", {
    obj_list <- make_seurat_list("Sample_A")
    sheet    <- make_sample_sheet("Sample_A", Treatment = "Control")

    res <- includeMetadata(obj_list, sheet)

    expect_named(res, "Sample_A")
    expect_equal(unique(res$Sample_A$Treatment), "Control")
})

test_that("overwrite = TRUE is harmless when nothing collides", {
    obj_list <- make_seurat_list(c("Sample_A", "Sample_B"))
    sheet    <- make_sample_sheet(c("Sample_A", "Sample_B"))

    res <- includeMetadata(obj_list, sheet, overwrite = TRUE)

    expect_equal(unique(res$Sample_A$Treatment), "Control")
})

test_that("integer columns stay integer", {
    obj_list <- make_seurat_list(c("Sample_A", "Sample_B"))
    sheet    <- make_sample_sheet(c("Sample_A", "Sample_B"),
                                  Passage = c(3L, 5L))

    res <- includeMetadata(obj_list, sheet)

    expect_type(res$Sample_A$Passage, "integer")
    expect_equal(unique(res$Sample_B$Passage), 5L)
})

test_that("columns are added alongside, not instead of, Seurat's own", {
    obj_list <- make_seurat_list(c("Sample_A", "Sample_B"))
    before   <- colnames(obj_list$Sample_A[[]])

    res <- includeMetadata(obj_list, make_sample_sheet(c("Sample_A", "Sample_B")))

    expect_true(all(before %in% colnames(res$Sample_A[[]])))
    expect_equal(unique(res$Sample_A$Ident), "Sample_A")
    expect_equal(colnames(res$Sample_A[[]]), c(before, "Treatment"))
})

test_that("annotating twice with different columns accumulates them", {
    obj_list <- make_seurat_list(c("Sample_A", "Sample_B"))

    res <- includeMetadata(obj_list, make_sample_sheet(c("Sample_A", "Sample_B")))
    res <- includeMetadata(res, make_sample_sheet(c("Sample_A", "Sample_B"),
                                                  Timepoint = c("D0", "D21")))

    expect_true(all(c("Treatment", "Timepoint") %in% colnames(res$Sample_A[[]])))
    expect_equal(unique(res$Sample_B$Timepoint), "D21")
})

# ---- integration: the reason completeness is enforced --------------------

test_that("annotated objects merge into one object with no NA cells", {
    root <- local_root()
    local_mock_read10x()
    outs <- make_multi_run(root, "Pool1", c("Sample_A", "Sample_B"))

    # one manifest row, two demuxed samples, two different treatments --
    # the case that cannot be expressed in the manifest itself
    seurat_list <- seuratFrom10x(path = outs, multi = TRUE,
                                 min.cells = 0, min.features = 0)
    sheet <- data.frame(Sample    = c("Sample_A", "Sample_B"),
                        Treatment = c("Control", "Bleomycin"))

    annotated <- includeMetadata(seurat_list, sheet)
    expect_output(merged <- mergeSeurat(annotated))

    expect_true("Treatment" %in% colnames(merged[[]]))
    expect_false(anyNA(merged$Treatment))
    expect_setequal(unique(merged$Treatment), c("Control", "Bleomycin"))
    expect_equal(unname(table(merged$Ident, merged$Treatment)["Sample_A", "Control"]), 21L)
})

test_that("partial annotation is what produces NA cells after a merge", {
    # Pins the Seurat behaviour includeMetadata() exists to prevent: merge()
    # accepts a ragged set of metadata columns without complaint.
    annotated   <- make_seurat("Sample_A", n_cells = 10)
    unannotated <- make_seurat("Sample_B", n_cells = 10)
    annotated$Treatment <- "Control"

    expect_output(merged <- mergeSeurat(list(Sample_A = annotated,
                                             Sample_B = unannotated)))

    expect_true("Treatment" %in% colnames(merged[[]]))
    expect_equal(sum(is.na(merged$Treatment)), 10L)
})
