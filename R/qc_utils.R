#== Internal helpers shared by the QC functions.
#==
#== diagnoseQC() is the only caller today, but the organism -> gene pattern
#== lookup and the "is this a usable Seurat object" check are both things any
#== future QC entry point will need, so they live here rather than inside the
#== exported function. Same split as data_loading_utils.R: the exported
#== function owns its argument surface, the workers here own the shared
#== machinery.

#== Regexes matching mitochondrial and ribosomal protein genes, keyed by the
#== three-letter organism code. These match gene *symbols*, and symbol case is
#== the one thing that differs between species: human symbols are upper case
#== throughout (MT-CO1, RPS6), mouse and zebrafish are title case with a lower
#== case prefix (mt-Co1, Rps6), and rat keeps the capitalised Mt- while its
#== ribosomal genes follow the mouse convention. The mouse mitochondrial
#== pattern is deliberately "^mt" without the hyphen -- some GTF versions emit
#== mt-Co1 and others mtCo1.
.organism_patterns <- list(
    hsa = list(mito = "^MT-",  ribo = "^RP[SL]"),
    mmu = list(mito = "^mt",   ribo = "^Rp[sl]"),
    rno = list(mito = "^Mt-",  ribo = "^Rp[sl]"),
    dre = list(mito = "^mt-",  ribo = "^rp[sl]")
)

#' Resolve an organism code to one of the supported set
#'
#' [match.arg()]'s own message is unreadable when the input is far off, and the
#' common failure -- a Latin name, a common name, or an Ensembl species code --
#' deserves a message naming what is actually accepted.
#'
#' @param organism The caller's argument, unmodified. A unique prefix is
#'   accepted, as is the full default vector.
#' @param call Environment to blame in the error, so the message points at the
#'   exported function rather than at this helper.
#'
#' @return A single supported organism code, or aborts.
#'
#' @noRd
resolveOrganism <- function(organism, call = parent.frame()) {
    supported <- names(.organism_patterns)

    if (!is.character(organism) || length(organism) == 0 || anyNA(organism)) {
        cli::cli_abort(
            "{.arg organism} must be a single organism code: {.val {supported}}.",
            call = call
        )
    }
    if (length(organism) > 1 && !identical(organism, supported)) {
        cli::cli_abort(
            "{.arg organism} must be a single organism code, not {length(organism)}.",
            call = call
        )
    }

    matched <- tryCatch(match.arg(organism, supported), error = function(e) NULL)
    if (is.null(matched)) {
        cli::cli_abort(c(
            "{.arg organism} {.val {organism}} is not supported.",
            "i" = "Supported organisms: {.val {supported}}.",
            "i" = "These are gene {.emph symbol} patterns, so an organism whose
                   symbols follow the same capitalisation convention as one of
                   these can borrow its code."
        ), call = call)
    }

    matched
}

#' Look up the mitochondrial and ribosomal gene patterns for an organism
#'
#' @param organism Three-letter organism code; one of `"hsa"` (human), `"mmu"`
#'   (mouse), `"rno"` (rat) or `"dre"` (zebrafish). Resolved by
#'   `resolveOrganism()`, so a unique prefix is accepted.
#' @param call Environment to blame in the error.
#'
#' @return A list of two regexes, `mito` and `ribo`, or aborts if `organism`
#'   is not one of the supported codes.
#'
#' @noRd
genePatterns <- function(organism = c("hsa", "mmu", "rno", "dre"),
                         call     = parent.frame()) {
    .organism_patterns[[resolveOrganism(organism, call = call)]]
}

#' Check that an object is a usable Seurat object for QC
#'
#' Structural validity only -- whether the metrics can be calculated at all.
#' Whether the metrics look *good* is what the QC functions are for.
#'
#' @param obj The caller's object, unmodified.
#' @param arg_name Name the caller passed it under, for the message.
#' @param fn Name of the exported function doing the checking, so the
#'   `lapply()` hint suggests the function the user actually called.
#' @param call Environment to blame in the error.
#'
#' @return `invisible(TRUE)`, or aborts.
#'
#' @importFrom methods validObject
#' @noRd
validateQCObject <- function(obj, arg_name, fn = "diagnoseQC",
                             call = parent.frame()) {
    #== The most likely mistake is handing over the named list seuratFrom10x()
    #== returns, so that case gets its own message rather than a bare class
    #== mismatch. data.frames are lists too, hence the explicit exclusion.
    if (is.list(obj) && !is.data.frame(obj) && !inherits(obj, "Seurat")) {
        cli::cli_abort(c(
            "{.arg {arg_name}} is a list, not a single Seurat object.",
            "i" = "To handle every sample in a list, either merge it first
                   with {.fn mergeSeurat} or map over it:
                   {.code lapply({arg_name}, {fn}, ...)}."
        ), call = call)
    }
    if (!inherits(obj, "Seurat")) {
        cli::cli_abort(
            "{.arg {arg_name}} must be a Seurat object, not {.cls {class(obj)[1]}}.",
            call = call
        )
    }

    valid <- tryCatch({
        methods::validObject(obj)
        TRUE
    }, error = function(e) FALSE)

    if (!valid) {
        cli::cli_abort("{.arg {arg_name}} is not a structurally valid Seurat object.",
                       call = call)
    }
    if (ncol(obj) == 0) {
        cli::cli_abort("{.arg {arg_name}} has zero cells -- nothing to diagnose.",
                       call = call)
    }
    if (nrow(obj) == 0) {
        cli::cli_abort("{.arg {arg_name}} has zero features -- nothing to diagnose.",
                       call = call)
    }

    invisible(TRUE)
}
