#' Add Per-Cell QC Metrics to a Seurat Object
#'
#' Calculates the three per-cell metrics that single-cell filtering decisions
#' are usually made from -- mitochondrial content, ribosomal content and
#' library complexity -- and writes them into the object's metadata as
#' `pctMito`, `pctRibo` and `log10GenesPerUMI`. Nothing is filtered and no
#' cells are dropped; the metrics are added so that thresholds can be chosen
#' by looking at the distributions first.
#'
#' The two percentage metrics are computed by
#' [Seurat::PercentageFeatureSet()] from organism-specific gene *symbol*
#' patterns, so the object's rownames must be symbols. A matrix still keyed by
#' Ensembl IDs will match nothing and produce columns of zeroes; that case is
#' warned about rather than silently returned, and the warning names the
#' pattern that failed so it is obvious what went wrong.
#'
#' @section The metrics:
#'
#' \describe{
#'   \item{`pctMito`}{Percentage of a cell's UMIs coming from mitochondrial
#'     genes. High values mark cells that lost cytoplasmic RNA through a
#'     ruptured membrane -- dying or stressed cells -- and are the usual first
#'     filter. What counts as high is tissue-dependent: 5-10% is a common
#'     starting point, but cardiac, muscle and kidney samples are legitimately
#'     mitochondria-rich and a fixed cutoff will throw away real cells.}
#'   \item{`pctRibo`}{Percentage coming from ribosomal protein genes. Less
#'     often used as a hard filter and more often as a covariate: it varies
#'     strongly by cell type, so a sharp shift between samples that should be
#'     comparable points at a technical difference rather than a biological
#'     one.}
#'   \item{`log10GenesPerUMI`}{Library complexity --
#'     `log10(nFeature) / log10(nCount)`. Falls between 0 and 1 for real
#'     cells, and low values mean a library dominated by few transcripts:
#'     ambient RNA, degraded cells, or red blood cells, which are genuinely
#'     low-complexity. Around 0.8 is a conventional floor.}
#' }
#'
#' The count columns the complexity metric is derived from
#' (`nCount_<assay>` / `nFeature_<assay>`) are read from the object's default
#' assay, so an object whose default assay is not `RNA` is handled without
#' any extra argument. The three metrics themselves are always named the same
#' way, because that is what the plotting functions look for.
#'
#' @section Gene patterns by organism:
#'
#' \describe{
#'   \item{`"hsa"` (human)}{mito `^MT-`, ribo `^RP[SL]`}
#'   \item{`"mmu"` (mouse)}{mito `^mt`, ribo `^Rp[sl]`}
#'   \item{`"rno"` (rat)}{mito `^Mt-`, ribo `^Rp[sl]`}
#'   \item{`"dre"` (zebrafish)}{mito `^mt-`, ribo `^rp[sl]`}
#' }
#'
#' These differ only in the capitalisation convention each species' gene
#' symbols follow, so an organism not listed here can borrow the code of one
#' whose convention it shares.
#'
#' @param obj A single [Seurat][Seurat::CreateSeuratObject] object, typically
#'   the merged output of [mergeSeurat()]. A named *list* of objects is not
#'   accepted -- merge it first, or `lapply()` over it -- so that the metrics
#'   are always calculated on the same object the QC plots are drawn from.
#' @param organism Three-letter code selecting the gene patterns: `"hsa"`
#'   (human, the default), `"mmu"` (mouse), `"rno"` (rat) or `"dre"`
#'   (zebrafish). See *Gene patterns by organism*.
#'
#' @return The input object with three numeric metadata columns added:
#'   `pctMito`, `pctRibo` and `log10GenesPerUMI`. Existing columns of the same
#'   name are overwritten, so the function is safe to re-run -- with a
#'   different `organism`, for instance. Every other slot is untouched, and
#'   the cell count is unchanged.
#'
#' @examples
#' \dontrun{
#' seurat_list <- seuratFrom10x(path = "/path/to/Pool1/outs", multi = TRUE)
#' seurat      <- mergeSeurat(seurat_list)
#'
#' seurat <- diagnoseQC(seurat, organism = "mmu")
#' head(seurat@meta.data)
#'
#' # or metric-by-sample, without merging
#' seurat_list <- lapply(seurat_list, diagnoseQC, organism = "mmu")
#' }
#'
#' @seealso [mergeSeurat()], which produces the object this expects;
#'   [plotBarcodeRanks()] for the other side of QC -- cell calling, from the
#'   raw matrix; [Seurat::PercentageFeatureSet()] for the percentage
#'   calculation itself.
#'
#' @export
diagnoseQC <- function(obj, organism = c("hsa", "mmu", "rno", "dre")) {
    arg_name <- deparse(substitute(obj))
    validateQCObject(obj, arg_name)

    #== Resolve before logging: the default is a 4-element vector, and an
    #== unresolved `organism` would report all four as the one in use.
    organism <- resolveOrganism(organism)
    patterns <- genePatterns(organism)

    #== Read the count columns off the default assay rather than hardcoding
    #== "RNA", so a re-keyed object still works without an extra argument.
    assay     <- Seurat::DefaultAssay(obj)
    count_col <- paste0("nCount_", assay)
    feat_col  <- paste0("nFeature_", assay)
    n_cells   <- ncol(obj)

    absent <- setdiff(c(count_col, feat_col), colnames(obj[[]]))
    if (length(absent) > 0) {
        cli::cli_abort(c(
            "{.arg {arg_name}} is missing the count {cli::qty(absent)}column{?s} {.val {absent}}.",
            "i" = "These are written by {.fn Seurat::CreateSeuratObject} for the
                   default assay, currently {.val {assay}}.",
            "i" = "Set {.fn Seurat::DefaultAssay} to the assay the counts came
                   from before diagnosing."
        ))
    }

    cli::cli_alert_info(
        "Diagnosing {n_cells} cell{?s} in assay {.val {assay}} using {.val {organism}} gene patterns."
    )

    #== Count matches before calculating, so a pattern matching nothing is
    #== reported as the cause rather than found later as a column of zeroes.
    metrics <- c(mito = "pctMito", ribo = "pctRibo")
    for (k in names(metrics)) {
        pattern  <- patterns[[k]]
        n        <- length(grep(pattern, rownames(obj)))
        col      <- metrics[[k]]
        peek     <- paste0("head(rownames(", arg_name, "))")

        if (n == 0) {
            cli::cli_warn(c(
                "!" = "{.field {k}} pattern {.val {pattern}} matched 0 features
                       for organism {.val {organism}}.",
                "i" = "Check that rownames are gene symbols and not Ensembl IDs
                       -- {.code {peek}} will show which.",
                "i" = "{.field {col}} will be 0 for every cell."
            ))
        } else {
            cli::cli_alert_info(
                "{.field {k}}: {n} feature{?s} matched {.val {pattern}}."
            )
        }
    }

    #== Percentage of each cell's counts coming from each feature set
    obj[["pctMito"]] <- Seurat::PercentageFeatureSet(obj, pattern = patterns$mito)
    obj[["pctRibo"]] <- Seurat::PercentageFeatureSet(obj, pattern = patterns$ribo)

    #== Library complexity: log10 genes per UMI
    complexity <- log10(obj[[feat_col]][[1]]) / log10(obj[[count_col]][[1]])
    obj[["log10GenesPerUMI"]] <- complexity

    #== A cell with a total count of 1 divides by log10(1) == 0. Rare, but it
    #== propagates into every downstream summary as NaN/Inf, so flag it once.
    finite <- is.finite(complexity)
    n_bad  <- sum(!finite)
    if (n_bad > 0) {
        cli::cli_warn(c(
            "!" = "{n_bad} cell{?s} {?has/have} a non-finite
                   {.field log10GenesPerUMI}.",
            "i" = "Caused by cells with a total count of 0 or 1, where the
                   denominator {.code log10(nCount)} is zero or undefined.",
            "i" = "They are left as-is; drop them with {.code is.finite()}
                   before summarising."
        ))
    }

    #== Summarise what was added, so an obviously wrong result (all-zero
    #== percentages, a complexity median far below 0.8) is visible immediately
    #== rather than only once the object is plotted.
    med_mito <- round(stats::median(obj$pctMito), 2)
    max_mito <- round(max(obj$pctMito), 2)
    med_ribo <- round(stats::median(obj$pctRibo), 2)
    max_ribo <- round(max(obj$pctRibo), 2)
    med_cplx <- if (any(finite)) round(stats::median(complexity[finite]), 3) else NA

    cli::cli_alert_success("Added 3 QC metrics to {n_cells} cell{?s}.")
    cli::cli_bullets(c(
        "*" = "{.field pctMito}: median {med_mito}%, max {max_mito}%",
        "*" = "{.field pctRibo}: median {med_ribo}%, max {max_ribo}%",
        "*" = "{.field log10GenesPerUMI}: median {med_cplx}"
    ))

    obj
}
