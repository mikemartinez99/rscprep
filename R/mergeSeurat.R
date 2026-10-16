#' Merge a Named List of Seurat Objects into One
#'
#' Combines a named list of per-sample
#' [Seurat][Seurat::CreateSeuratObject] objects -- typically the output of
#' [seuratFrom10x()] -- into a single object via [merge()], prefixing every
#' cell barcode with its sample name so barcodes stay unique and traceable
#' back to their sample of origin.
#'
#' The list names are load-bearing: they are passed to `merge()` as
#' `add.cell.id`, so a cell `AAACCCAAGA-1` in list entry `Sample_1` becomes
#' `Sample_1_AAACCCAAGA-1` in the merged object. Names must therefore be
#' present, non-blank and unique.
#'
#' The whole list is validated up front, so a malformed batch fails with a
#' specific message instead of partway through the merge. In order, `objList`
#' must:
#'
#' * be non-empty (a single-element list short-circuits -- see *Value*);
#' * contain only structurally valid Seurat objects;
#' * be fully named, with no `NA`, blank or duplicated names;
#' * contain no zero-cell objects;
#' * carry an `Ident` metadata column, which [seuratFrom10x()] sets to the
#'   sample name;
#' * agree on [Seurat::DefaultAssay()] across every element.
#'
#' Two further checks run *after* the merge, guarding against a merge that
#' silently loses data: every sample must still be represented in the merged
#' `Ident` metadata column, and the merged cell count must equal the sum of
#' the inputs'.
#'
#' @section Seurat 4 and Seurat 5:
#' The installed version is detected automatically. Under Seurat 4, `merge()`
#' binds every sample's counts into one matrix. Under Seurat 5 it keeps one
#' layer per sample (`counts.Sample_A`, `counts.Sample_B`, ...), which is
#' what Seurat 5's `IntegrateLayers()` wants but not what most other code
#' expects. By default those layers are joined, so the merged object has the
#' same shape whichever version made it. Set `joinLayers = FALSE` to keep them
#' split; under Seurat 4 the argument has no effect.
#'
#' @param objList Named list of [Seurat][Seurat::CreateSeuratObject] objects,
#'   one per sample, as returned by [seuratFrom10x()]. Names are used as cell
#'   barcode prefixes and must be unique.
#' @param joinLayers Under Seurat 5, join the per-sample layers `merge()`
#'   creates back into one per kind of data. Defaults to `TRUE`. Ignored under
#'   Seurat 4, which never splits them. See *Seurat 4 and Seurat 5*.
#'
#' @return A single merged [Seurat][Seurat::CreateSeuratObject] object,
#'   returned **invisibly** (the object is printed as a side effect). Cell
#'   barcodes are prefixed with their sample name; the `Ident` column is
#'   carried through from the inputs. If `objList` has exactly one element,
#'   that element is returned unmerged and unprefixed, with a message.
#'
#' @examples
#' \dontrun{
#' seurat_list <- seuratFrom10x(
#'   path        = "/path/to/Pool1/outs",
#'   multi = TRUE
#' )
#'
#' merged <- mergeSeurat(seurat_list)
#'
#' # barcodes are namespaced by sample
#' head(colnames(merged))
#'
#' # and Ident still identifies the sample of origin
#' table(merged$Ident)
#' }
#'
#' @seealso [seuratFrom10x()], which produces the named list this expects.
#'
#' @importFrom Seurat DefaultAssay
#' @importFrom methods validObject
#' @export
mergeSeurat <- function(objList, joinLayers = TRUE) {
    arg_name <- deparse(substitute(objList))

    #== Guard against empty or trivial input
    if (length(objList) == 0) {
        cli::cli_abort("'{arg_name}' is empty -- nothing to merge.")
    }
    if (length(objList) == 1) {
        cli::cli_alert_info("'{arg_name}' has only one object; returning it unmerged.")
        return(objList[[1]])
    }

    #== Validate each element is a real, structurally valid Seurat object
    is_valid <- vapply(objList, function(x) {
        inherits(x, "Seurat") && tryCatch({
            methods::validObject(x)
            TRUE
        }, error = function(e) FALSE)
    }, logical(1))

    if (!all(is_valid)) {
        bad <- if (!is.null(names(objList))) names(objList)[!is_valid] else which(!is_valid)
        cli::cli_abort("The following element(s) of '{arg_name}' are not valid Seurat objects: {bad}")
    }

    samples <- names(objList)
    n_objects <- length(objList)

    #== Names must exist, be non-blank, and be unique -- add.cell.id relies on this
    if (is.null(samples) || anyNA(samples) || any(samples == "")) {
        cli::cli_abort("'{arg_name}' must be fully named -- found missing or blank name(s).")
    }
    if (anyDuplicated(samples)) {
        dupes <- unique(samples[duplicated(samples)])
        cli::cli_abort("'{arg_name}' contains duplicate names: {dupes}")
    }

    #== No zero-cell objects
    n_cells <- vapply(objList, ncol, numeric(1))   # double under SeuratObject 5
    if (any(n_cells == 0)) {
        empty_samples <- samples[n_cells == 0]
        cli::cli_abort("The following sample(s) have zero cells: {empty_samples}")
    }

    #== Every object must carry Ident -- the post-merge check reads it. Checked
    #== here because a missing column otherwise surfaces as SeuratObject's own
    #== error from `merged_obj$Ident`, worded differently by Seurat 4 and 5
    #== and naming neither the argument nor the sample.
    no_ident <- samples[!vapply(objList, function(x) "Ident" %in% colnames(x@meta.data),
                                logical(1))]
    if (length(no_ident) > 0) {
        cli::cli_abort(c(
            "The following element(s) of '{arg_name}' have no {.field Ident} column: {no_ident}",
            "i" = "{.fn seuratFrom10x} sets it to the sample name; set it the
                   same way on objects built another way."
        ))
    }

    #== All objects must share the same default assay
    default_assays <- vapply(objList, Seurat::DefaultAssay, character(1))
    if (length(unique(default_assays)) > 1) {
        cli::cli_abort(
            "'{arg_name}' contains objects with mismatched default assays: {unique(default_assays)}"
        )
    }

    cli::cli_alert_info("Merging {n_objects} seurat objects...")
    merged_obj <- merge(x = objList[[1]],
                        y = objList[-1],
                        add.cell.id = samples)

    #== Every sample must survive the merge
    if (n_objects != length(unique(merged_obj$Ident))) {
        missing <- setdiff(samples, unique(merged_obj$Ident))
        cli::cli_abort("The following element(s) of '{arg_name}' were not merged correctly: {missing}")
    }

    #== No cells silently dropped or duplicated during the merge
    if (sum(n_cells) != ncol(merged_obj)) {
        cli::cli_abort(
            "Cell count mismatch after merging '{arg_name}': expected {sum(n_cells)}, got {ncol(merged_obj)}."
        )
    }

    #== Seurat 5 only: merge() left one layer per sample; join them so the
    #== result matches what Seurat 4 returns. A no-op under Seurat 4.
    if (isTRUE(joinLayers)) merged_obj <- joinAssayLayers(merged_obj)

    cli::cli_alert_success("Successfully merged {n_objects} seurat objects.")
    print(merged_obj)
    return(invisible(merged_obj))
}
