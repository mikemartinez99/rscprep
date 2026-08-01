#' Attach Sample-Level Metadata to a List of Seurat Objects
#'
#' Joins a sample-level table -- treatment group, genotype, timepoint, and so
#' on -- onto a named list of [Seurat][Seurat::CreateSeuratObject] objects,
#' typically the output of [seuratFrom10x()]. Each row of `metadata`
#' describes one sample; every other column becomes a metadata column on that
#' sample's object, broadcast across all of its cells.
#'
#' This is deliberately separate from the manifest [seuratFrom10x()] accepts.
#' That manifest is a *run*-level table -- one row per Cell Ranger `outs/`
#' directory -- and a multiplexed run expands into many samples, so a single
#' manifest row cannot carry per-sample assignments. `metadata` is a
#' *sample*-level table instead, keyed on the names of `objList`. That works
#' identically for multiplexed and non-multiplexed input, because by the time
#' the list exists the distinction is already gone.
#'
#' Samples are identified by a `Sample` column if `metadata` has one, and by
#' its rownames otherwise. A `Sample` column takes precedence, and is used as
#' the key rather than attached as metadata.
#'
#' Every column must be populated for *every* object in `objList`: a sample
#' with no row in `metadata`, or an `NA` or blank value in any column, is an
#' error. This is what keeps the result safe to [mergeSeurat()]. `merge()`
#' will happily combine an object carrying a `Treatment` column with one that
#' lacks it, silently filling the gap with `NA` for every cell of the
#' un-annotated sample; requiring completeness up front makes that state
#' unreachable.
#'
#' Rows matching no object are ignored with a warning, so one master sheet can
#' cover a whole experiment while you work on a subset of it. A mistyped
#' sample name is still caught, since it also leaves some object without a row.
#'
#' @param objList Named list of [Seurat][Seurat::CreateSeuratObject] objects,
#'   one per sample, as returned by [seuratFrom10x()].
#' @param metadata Data frame with one row per sample. Sample names come from
#'   a `Sample` column if present, otherwise from the rownames, and must match
#'   `names(objList)`. Every remaining column is attached.
#' @param overwrite Logical. Whether to replace metadata columns that already
#'   exist on the objects. Defaults to `FALSE`, which errors instead. `Ident`
#'   is reserved and cannot be set either way -- it identifies the sample of
#'   origin and is what [mergeSeurat()] verifies after merging.
#'
#' @return `objList` with the new columns attached, in the same order and with
#'   the same names. Column types are preserved, including factor levels --
#'   though note that `Seurat::merge()` flattens factors to character, so a
#'   factor column survives on the list but not through [mergeSeurat()].
#'
#' @examples
#' \dontrun{
#' seurat_list <- seuratFrom10x(
#'   path        = "/path/to/Pool1/outs",
#'   multiplexed = TRUE
#' )
#'
#' # one row per demuxed sample, however many the run produced
#' sheet <- data.frame(
#'   Sample    = c("Sample_A", "Sample_B"),
#'   Treatment = c("Control", "Bleomycin"),
#'   Timepoint = c("D0", "D21")
#' )
#'
#' annotated <- includeMetadata(seurat_list, sheet)
#' merged    <- mergeSeurat(annotated)
#'
#' table(merged$Ident, merged$Treatment)
#' }
#'
#' @seealso [seuratFrom10x()], which produces the named list this expects, and
#'   [mergeSeurat()], which consumes the annotated result.
#'
#' @export
includeMetadata <- function(objList, metadata, overwrite = FALSE) {

    #== objList must be a named list of Seurat objects, not a single object
    if (inherits(objList, "Seurat")) {
        cli::cli_abort(c(
            "{.arg objList} must be a named {.cls list} of Seurat objects, not a single object.",
            "i" = "Metadata is attached per sample, before merging."
        ))
    }
    if (!is.list(objList)) {
        cli::cli_abort(
            "{.arg objList} must be a named {.cls list} of Seurat objects, not {.cls {class(objList)[1]}}."
        )
    }
    if (length(objList) == 0) {
        cli::cli_abort("{.arg objList} is empty -- nothing to annotate.")
    }

    samples <- names(objList)
    if (is.null(samples) || anyNA(samples) || any(samples == "")) {
        cli::cli_abort("{.arg objList} must be fully named -- found missing or blank names.")
    }
    if (anyDuplicated(samples)) {
        dupes <- unique(samples[duplicated(samples)])
        cli::cli_abort("{.arg objList} contains duplicate names: {.val {dupes}}")
    }

    not_seurat <- !vapply(objList, inherits, logical(1), "Seurat")
    if (any(not_seurat)) {
        bad <- samples[not_seurat]
        cli::cli_abort(
            "{length(bad)} element{?s} of {.arg objList} {?is/are} not a {.cls Seurat} object: {.val {bad}}"
        )
    }

    #== metadata must be a rectangular table with rows
    if (!is.data.frame(metadata)) {
        cli::cli_abort(
            "{.arg metadata} must be a {.cls data.frame}, not {.cls {class(metadata)[1]}}."
        )
    }
    metadata <- as.data.frame(metadata)   # tibbles, data.tables, ...
    if (nrow(metadata) == 0) {
        cli::cli_abort("{.arg metadata} has no rows.")
    }

    #== Sample column wins; rownames are the fallback key
    if ("Sample" %in% colnames(metadata)) {
        keys      <- as.character(metadata[["Sample"]])
        meta_cols <- setdiff(colnames(metadata), "Sample")
    } else {
        # negative row.names info == automatic 1:n rownames, which would
        # silently be matched against sample names as "1", "2", "3"
        if (.row_names_info(metadata) < 0) {
            cli::cli_abort(c(
                "Can't tell which sample each row of {.arg metadata} describes.",
                "i" = "Add a {.field Sample} column, or set the rownames to the sample names.",
                "i" = "Either way they must match {.code names(objList)}: {.val {samples}}"
            ))
        }
        keys      <- rownames(metadata)
        meta_cols <- colnames(metadata)
    }

    if (length(meta_cols) == 0) {
        cli::cli_abort(c(
            "{.arg metadata} has no metadata columns -- only the {.field Sample} key.",
            "i" = "Add at least one column to attach."
        ))
    }

    # Seurat silently runs make.names() on metadata column names, so
    # "Cell Type" would land as `Cell.Type` and no longer match the sheet
    # header the user is reading from. Refuse rather than rename behind them.
    invalid_cols <- meta_cols[meta_cols != make.names(meta_cols)]
    if (length(invalid_cols) > 0) {
        cli::cli_abort(c(
            "{length(invalid_cols)} column name{?s} in {.arg metadata} {?is/are} not syntactically valid: {.val {invalid_cols}}",
            "i" = "Seurat would silently rename {.val {invalid_cols}} to {.val {make.names(invalid_cols)}}.",
            "i" = "Rename them in {.arg metadata} first."
        ))
    }
    if (anyNA(keys) || any(trimws(keys) %in% "")) {
        cli::cli_abort("{.arg metadata} contains missing or blank sample names.")
    }
    if (anyDuplicated(keys)) {
        dupes <- unique(keys[duplicated(keys)])
        cli::cli_abort("{.arg metadata} contains duplicate sample names: {.val {dupes}}")
    }

    #== Every sample needs a row -- this is half of the completeness contract
    missing_samples <- setdiff(samples, keys)
    if (length(missing_samples) > 0) {
        cli::cli_abort(c(
            "Every object in {.arg objList} needs a row in {.arg metadata}.",
            "x" = "No row for: {.val {missing_samples}}",
            "i" = "Metadata must cover every sample so it stays consistent after {.fun mergeSeurat}."
        ))
    }

    #== Rows for samples we do not have are harmless -- say so and move on
    extra_rows <- setdiff(keys, samples)
    if (length(extra_rows) > 0) {
        cli::cli_warn(
            "Ignoring {length(extra_rows)} row{?s} of {.arg metadata} matching no object in {.arg objList}: {.val {extra_rows}}"
        )
    }

    metadata <- metadata[match(samples, keys), meta_cols, drop = FALSE]
    rownames(metadata) <- samples

    #== ...and the other half: no gaps within the columns we do have
    incomplete <- lapply(meta_cols, function(col) samples[isBlank(metadata[[col]])])
    names(incomplete) <- meta_cols
    incomplete <- incomplete[lengths(incomplete) > 0]

    if (length(incomplete) > 0) {
        # report every offending column at once -- fixing a sheet one error
        # per run is miserable at 40 samples
        bullets <- vapply(names(incomplete), function(col) {
            cli::format_inline("{.field {col}}: missing for {.val {incomplete[[col]]}}")
        }, character(1), USE.NAMES = FALSE)
        # cli_abort() glue-interpolates each bullet, and these are already
        # formatted -- a brace in a sample name would be evaluated as an
        # expression and blow up the error path itself
        bullets <- escapeBraces(bullets)
        names(bullets) <- rep("x", length(bullets))

        cli::cli_abort(c(
            "Every column of {.arg metadata} must be populated for every sample.",
            bullets,
            "i" = "A column left blank for some samples becomes {.code NA} for their cells after {.fun mergeSeurat}."
        ))
    }

    #== Collisions with metadata the objects already carry
    if ("Ident" %in% meta_cols) {
        cli::cli_abort(c(
            "{.field Ident} is reserved and cannot be set from {.arg metadata}.",
            "i" = "It identifies the sample of origin, and is what {.fun mergeSeurat} verifies after merging."
        ))
    }

    existing <- unique(unlist(lapply(objList, function(obj) colnames(obj[[]]))))
    clashes  <- intersect(meta_cols, existing)
    if (length(clashes) > 0 && !overwrite) {
        cli::cli_abort(c(
            "{length(clashes)} column{?s} in {.arg metadata} {?is/are} already present in {.arg objList}: {.val {clashes}}",
            "i" = "Set {.code overwrite = TRUE} to replace them."
        ))
    }

    #== Broadcast each sample's row across all of its cells
    objList[] <- lapply(samples, function(sample) {
        obj <- objList[[sample]]
        for (col in meta_cols) {
            obj[[col]] <- rep(metadata[sample, col], ncol(obj))
        }
        obj
    })

    cli::cli_alert_success(
        "Added {length(meta_cols)} metadata column{?s} to {length(objList)} object{?s}: {.val {meta_cols}}"
    )

    objList
}

#' Escape braces in already-formatted cli text.
#'
#' [cli::cli_abort()] glue-interpolates every element it is handed, including
#' bullets built earlier with [cli::format_inline()]. Doubling the braces
#' makes glue emit them literally instead of evaluating what is between them.
#'
#' @noRd
escapeBraces <- function(x) {
    gsub("}", "}}", gsub("{", "{{", x, fixed = TRUE), fixed = TRUE)
}

#' Is a value missing, for metadata purposes?
#'
#' `NA` for any type, plus empty or whitespace-only strings -- a blank cell in
#' a CSV reads as `""`, not `NA`, and is just as unusable downstream.
#'
#' @noRd
isBlank <- function(x) {
    blank <- is.na(x)
    if (is.character(x) || is.factor(x)) {
        blank <- blank | trimws(as.character(x)) == ""
    }
    blank
}
