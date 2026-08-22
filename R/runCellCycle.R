#' Score Cells by Cell Cycle Phase
#'
#' Assigns each cell an S-phase score, a G2/M score and a phase call, using the
#' Tirosh et al. marker sets that ship with \pkg{Seurat} and
#' [Seurat::CellCycleScoring()]. The marker lists are human symbols; for any
#' other organism they are converted to that species' orthologs first, so the
#' same call works across species without the caller assembling gene lists.
#'
#' Cell cycle is a QC concern because it is one of the largest sources of
#' variation in a proliferating population, and it is *structured* variation --
#' cycling cells cluster with cycling cells regardless of type. The scores are
#' added here so that they can be looked at, and if necessary regressed out
#' later; nothing is regressed, scaled or clustered by this function.
#'
#' @section What gets added:
#'
#' [Seurat::CellCycleScoring()] writes three metadata columns:
#'
#' \describe{
#'   \item{`S.Score`, `G2M.Score`}{Mean expression of each marker set, centred
#'     against a background of control features drawn from comparable
#'     expression bins. Relative, not absolute: a score near zero means "no
#'     more than the background", not "no expression".}
#'   \item{`Phase`}{`"S"`, `"G2M"` or `"G1"`, from whichever score is higher,
#'     with `"G1"` when both are negative. `"G1"` therefore also absorbs G0 --
#'     the two are not separable from these markers.}
#' }
#'
#' `Idents()` is deliberately left alone (`set.ident = FALSE`), so the sample
#' identity [seuratFrom10x()] set survives.
#'
#' @section Normalisation:
#'
#' Scoring bins genes by average expression, which assumes the data have been
#' normalised. An object that still looks like raw counts is normalised here
#' with [Seurat::NormalizeData()] at its defaults, after a warning -- and the
#' returned object carries that normalised data, since scoring against one
#' version of the data and handing back another would be worse than either.
#'
#' Normalise beforehand if you want anything other than the defaults; the check
#' is a heuristic (log-normalised values are essentially never whole numbers),
#' so an already-normalised object is left exactly as it was.
#'
#' @section Non-human organisms:
#'
#' Markers are converted with [orthogene::convert_orthologs()], which is only
#' needed -- and only loaded -- when `organism` is not `"hsa"`. Genes without a
#' one-to-one ortholog are dropped, so the converted sets are smaller than the
#' human ones; how many survived is reported.
#'
#' The default conversion `method`, `"gprofiler"`, queries a web service and so
#' needs network access. On a machine without it, pass
#' `method = "homologene"`, which resolves locally from a bundled mapping at
#' the cost of some coverage.
#'
#' @param obj A single [Seurat][Seurat::CreateSeuratObject] object, normalised.
#' @param organism Three-letter code selecting the marker set: `"hsa"` (human,
#'   the default) uses [Seurat::cc.genes] directly; `"mmu"` (mouse), `"rno"`
#'   (rat) and `"dre"` (zebrafish) convert it to orthologs first.
#' @param method Ortholog mapping backend passed to
#'   [orthogene::convert_orthologs()]. `"gprofiler"` (the default) is the most
#'   complete but needs network access; `"homologene"` and `"babelgene"`
#'   resolve locally. Ignored when `organism = "hsa"`.
#'
#' @return The input object with `S.Score`, `G2M.Score` and `Phase` added to
#'   its metadata. Cell count and identities are unchanged. The data layer is
#'   unchanged too, unless the object arrived unnormalised -- see
#'   *Normalisation*.
#'
#' @examples
#' \dontrun{
#' seurat <- mergeSeurat(seuratFrom10x(path = "/path/to/Pool1/outs", multi = TRUE))
#' seurat <- Seurat::NormalizeData(seurat)
#'
#' seurat <- runCellCycle(seurat, organism = "mmu")
#'
#' table(seurat$Phase)
#'
#' # is one sample cycling harder than the others?
#' table(seurat$Ident, seurat$Phase)
#'
#' # on a cluster with no outbound network
#' seurat <- runCellCycle(seurat, organism = "mmu", method = "homologene")
#'
#' # the scores are ordinary metadata, so they plot like any other metric
#' plotScatter(seurat, "S.Score", "G2M.Score", group_by = "Phase")
#' }
#'
#' @seealso [Seurat::CellCycleScoring()] for the scoring itself, [diagnoseQC()]
#'   for the other per-cell QC metrics.
#'
#' @export
runCellCycle <- function(obj,
                         organism = c("hsa", "mmu", "rno", "dre"),
                         method   = c("gprofiler", "homologene", "babelgene")) {

    arg_name <- deparse(substitute(obj))
    validateQCObject(obj, arg_name, fn = "runCellCycle")

    #== Resolve before logging: the default is a 4-element vector, and an
    #== unresolved `organism` would report all four as the one in use.
    organism <- resolveOrganism(organism)
    method   <- match.arg(method)

    #== Scoring bins genes by average expression, which is meaningless on raw
    #== counts, so an unnormalised object is normalised here rather than
    #== scored as-is. It warns first: this changes the data layer of the object
    #== the caller gets back, which is not something to do quietly.
    if (isTRUE(looksUnnormalised(obj))) {
        cli::cli_warn(c(
            "!" = "{.arg {arg_name}} looks like it still holds raw counts.",
            "i" = "{.fn Seurat::CellCycleScoring} bins genes by average
                   expression, so {.fn Seurat::NormalizeData} is being run
                   first with its defaults.",
            "i" = "The returned object carries that normalised data. Normalise
                   it yourself beforehand if you want different parameters."
        ))
        obj <- tryCatch(
            Seurat::NormalizeData(obj, verbose = FALSE),
            error = function(e) {
                cli::cli_abort(c(
                    "{.fn Seurat::NormalizeData} failed on {.arg {arg_name}}.",
                    "x" = escapeBraces(conditionMessage(e)),
                    "i" = "Scoring needs normalised data; normalise the object
                           yourself and call {.fn runCellCycle} again."
                ))
            }
        )
    }

    features <- cellCycleFeatures(organism, method)
    s_genes   <- features$s
    g2m_genes <- features$g2m

    #== The markers have to be present in this object, not merely exist for the
    #== species -- the usual cause of an empty overlap is the wrong organism.
    #== Subsetting to what is present also stops CellCycleScoring() printing a
    #== list of every absent gene, and keeps the counts reported here equal to
    #== the counts actually scored.
    s_genes   <- s_genes[s_genes %in% rownames(obj)]
    g2m_genes <- g2m_genes[g2m_genes %in% rownames(obj)]
    n_s       <- length(s_genes)
    n_g2m     <- length(g2m_genes)

    if (n_s == 0 || n_g2m == 0) {
        cli::cli_abort(c(
            "None of the {.val {organism}} {if (n_s == 0) 'S' else 'G2/M'}
             markers are present in {.arg {arg_name}}.",
            "i" = "Matched {n_s} of {length(features$s)} S and {n_g2m} of
                   {length(features$g2m)} G2/M markers.",
            "i" = "Check that {.arg organism} matches the data and that
                   rownames are gene symbols, not Ensembl IDs."
        ))
    }

    cli::cli_alert_info(
        "Scoring {ncol(obj)} cell{?s} against {n_s} S and {n_g2m} G2/M \\
         {.val {organism}} markers."
    )
    if (n_s < 10 || n_g2m < 10) {
        cli::cli_warn(c(
            "!" = "Few markers matched: {n_s} S, {n_g2m} G2/M.",
            "i" = "Scores from a handful of genes are dominated by those genes
                   rather than by the phase they stand for."
        ))
    }

    #== CellCycleScoring() fails on small or low-complexity objects, with a
    #== message about bins or about sampling that says nothing about why.
    obj <- tryCatch(
        Seurat::CellCycleScoring(obj,
                                 s.features   = s_genes,
                                 g2m.features = g2m_genes,
                                 set.ident    = FALSE),
        error = function(e) {
            cli::cli_abort(c(
                "{.fn Seurat::CellCycleScoring} failed.",
                "x" = escapeBraces(conditionMessage(e)),
                "i" = "Complaints about bins or about sampling mean there were
                       too few genes or cells to build the control set that
                       scoring compares the markers against."
            ))
        }
    )

    phases <- table(factor(obj$Phase, levels = c("G1", "S", "G2M")))
    cli::cli_alert_success(
        "Added {.field S.Score}, {.field G2M.Score} and {.field Phase}."
    )
    cli::cli_bullets(c(
        "*" = "{.field G1}: {phases[['G1']]} cell{?s}",
        "*" = "{.field S}: {phases[['S']]} cell{?s}",
        "*" = "{.field G2M}: {phases[['G2M']]} cell{?s}"
    ))

    obj
}

#== orthogene's species names, keyed by the same three-letter codes the rest
#== of the package uses.
.cc_species <- c(hsa = "human", mmu = "mouse", rno = "rat", dre = "zebrafish")

#' The S and G2/M marker sets for one organism
#'
#' Human is the identity case: [Seurat::cc.genes] is already human symbols.
#' Everything else goes through the ortholog conversion.
#'
#' @param organism A resolved three-letter organism code.
#' @param method Backend passed to [orthogene::convert_orthologs()].
#' @param call Environment to blame in the error.
#'
#' @return A list of two character vectors, `s` and `g2m`.
#'
#' @noRd
cellCycleFeatures <- function(organism, method, call = parent.frame()) {
    s   <- Seurat::cc.genes$s.genes
    g2m <- Seurat::cc.genes$g2m.genes

    if (organism == "hsa") {
        return(list(s = s, g2m = g2m))
    }
    cellCycleOrthologs(s, g2m, organism, method, call = call)
}

#' Convert the human cell cycle markers to another species' orthologs
#'
#' Both marker sets are converted in a single [orthogene::convert_orthologs()]
#' call and split apart afterwards, rather than one call each: with the
#' default `gprofiler` backend each call is a network round trip, and the two
#' sets do not overlap, so one request does the work of two.
#'
#' Results are matched back by input symbol rather than by position, because
#' genes without a one-to-one ortholog are dropped by the conversion and the
#' output is neither the same length nor the same order as the input.
#'
#' @param s_features,g2m_features Human gene symbols for each phase.
#' @param organism A resolved three-letter organism code, not `"hsa"`.
#' @param method Backend passed to [orthogene::convert_orthologs()].
#' @param call Environment to blame in the error.
#'
#' @return A list of two character vectors, `s` and `g2m`, in the target
#'   species' symbols, with unmapped genes dropped.
#'
#' @noRd
cellCycleOrthologs <- function(s_features, g2m_features, organism, method,
                               call = parent.frame()) {
    if (!requireNamespace("orthogene", quietly = TRUE)) {
        cli::cli_abort(c(
            "The {.pkg orthogene} package is required to score a
             {.val {organism}} object.",
            "i" = "Install it with
                   {.code BiocManager::install(\"orthogene\")}, or score a
                   human object, which needs no conversion."
        ), call = call)
    }

    species <- .cc_species[[organism]]
    input   <- unique(c(s_features, g2m_features))

    cli::cli_alert_info(
        "Converting {length(input)} human cell cycle marker{?s} to \\
         {.val {species}} orthologs via {.val {method}}."
    )

    converted <- tryCatch(
        orthogene::convert_orthologs(
            gene_df        = input,
            gene_output    = "columns",
            input_species  = "human",
            output_species = species,
            method         = method,
            verbose        = FALSE
        ),
        error = function(e) {
            cli::cli_abort(c(
                "Ortholog conversion to {.val {species}} failed.",
                "x" = escapeBraces(conditionMessage(e)),
                "i" = "{.val {method}} needs network access if it is
                       {.val gprofiler}; {.val homologene} resolves locally."
            ), call = call)
        }
    )

    #== Genes without a one-to-one ortholog are dropped, so the result is
    #== neither the same length nor the same order as the input: match on the
    #== input symbol rather than trusting position.
    map <- stats::setNames(as.character(converted$ortholog_gene),
                           as.character(converted$input_gene))

    pick <- function(x) unname(map[x][!is.na(map[x])])
    out  <- list(s = pick(s_features), g2m = pick(g2m_features))

    if (length(out$s) == 0 || length(out$g2m) == 0) {
        cli::cli_abort(c(
            "No {.val {species}} orthologs were found for the
             {if (length(out$s) == 0) 'S' else 'G2/M'} markers.",
            "i" = "Converted {length(out$s)} of {length(s_features)} S and
                   {length(out$g2m)} of {length(g2m_features)} G2/M markers."
        ), call = call)
    }

    n_dropped <- length(s_features) + length(g2m_features) -
        length(out$s) - length(out$g2m)
    if (n_dropped > 0) {
        cli::cli_alert_info(
            "{n_dropped} marker{?s} had no one-to-one {.val {species}} ortholog \\
             and {cli::qty(n_dropped)}{?was/were} dropped."
        )
    }

    out
}

#' Does the default assay's data layer still look like raw counts?
#'
#' A heuristic, and deliberately a cheap one: log-normalised values are almost
#' never whole numbers, so all-integer data means normalisation has not run.
#' Only the first few thousand non-zero values are examined, and anything that
#' goes wrong reading the assay returns `NA` so the caller skips the warning
#' rather than failing on a Seurat version difference.
#'
#' @param obj A Seurat object.
#'
#' @return `TRUE`, `FALSE`, or `NA` when the check could not be made.
#'
#' @noRd
looksUnnormalised <- function(obj) {
    tryCatch({
        d <- Seurat::GetAssayData(obj, slot = "data")
        x <- if (methods::is(d, "sparseMatrix")) d@x else as.vector(d)
        x <- x[x != 0]

        if (length(x) == 0) return(NA)
        if (length(x) > 5000) x <- x[seq_len(5000)]

        all(x == round(x))
    }, error = function(e) NA)
}
