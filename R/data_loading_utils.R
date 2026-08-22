#== Internal workers shared by seuratFrom10x(), rawMatrixFrom10x() and
#== plotBarcodeRanks().
#==
#== Two things live here. First, the argument validation both entry points
#== share -- they take the same `manifest`/`path`/`multi` triple, so the
#== checks and their messages are defined once. Second, the actual file I/O.
#==
#== Every reader takes an `on_missing` callback rather than branching on a mode
#== flag: the default warns and returns NULL (manifest mode, where one absent
#== sample should not kill a batch), and `path` mode injects a handler that
#== aborts instead.
#==
#== All three entry points call prepareSource() rather than the individual
#== validators, so the deprecated `multiplexed` spelling is translated in one
#== place instead of three.

#' Validate the `manifest` / `path` / `multi` triple
#'
#' @param manifest,path,multi The caller's arguments, unmodified.
#' @param call Environment to blame in the error, so the message points at the
#'   exported function rather than at this helper.
#'
#' @return `invisible(TRUE)`, or aborts.
#'
#' @noRd
validateSourceArgs <- function(manifest, path, multi, call = parent.frame()) {
    if (is.null(manifest) && is.null(path)) {
        cli::cli_abort("Provide either {.arg manifest} or {.arg path}.", call = call)
    }
    if (!is.null(manifest) && !is.null(path)) {
        cli::cli_abort("Provide only one of {.arg manifest} or {.arg path}, not both.",
                       call = call)
    }
    if (!is.null(path) && is.null(multi)) {
        cli::cli_abort(
            "When using {.arg path}, you must also specify {.arg multi} ({.val {TRUE}}/{.val {FALSE}}).",
            call = call
        )
    }
    invisible(TRUE)
}

#' Translate the deprecated `multiplexed` spelling
#'
#' The flag was renamed because it selects a Cell Ranger *output layout*, not a
#' property of the experiment. `cellranger multi` writes `per_sample_outs/`
#' whether it demuxed twelve samples or one, so a singleplex Flex or 5' run
#' still needs `multi = TRUE` -- which reads as a contradiction when the
#' argument is called `multiplexed`.
#'
#' @param multi,multiplexed The new and old arguments, as supplied.
#' @param call Environment to blame in the error.
#'
#' @return The resolved value, possibly `NULL`.
#'
#' @noRd
resolveMultiArg <- function(multi, multiplexed, call = parent.frame()) {
    if (is.null(multiplexed)) {
        return(multi)
    }
    if (!is.null(multi)) {
        cli::cli_abort(
            c("Provide only one of {.arg multi} or {.arg multiplexed}, not both.",
              "i" = "{.arg multiplexed} is the old name for {.arg multi}."),
            call = call
        )
    }
    cli::cli_warn(c(
        "{.arg multiplexed} is deprecated; use {.arg multi} instead.",
        "i" = "It asks whether the run was processed with {.code cellranger multi}, which is also how singleplex Flex and 5' runs are processed."
    ))
    multiplexed
}

#' Translate a manifest that still uses the `multiplexed` column
#'
#' @param manifest The manifest, as supplied.
#' @param call Environment to blame in the error.
#'
#' @return The manifest, with a `multi` column.
#'
#' @noRd
migrateManifestCols <- function(manifest, call = parent.frame()) {
    #== Non-frames fall through untouched so validateManifest() reports the
    #== missing columns rather than this helper failing on colnames(NULL).
    cols <- colnames(manifest)
    if ("multi" %in% cols || !("multiplexed" %in% cols)) {
        return(manifest)
    }

    cli::cli_warn(c(
        "The {.field multiplexed} manifest column is deprecated; rename it to {.field multi}.",
        "i" = "It asks whether the run was processed with {.code cellranger multi}, which is also how singleplex Flex and 5' runs are processed."
    ))
    manifest$multi <- manifest$multiplexed
    manifest
}

#' Validate and normalise the source arguments every entry point shares
#'
#' One call instead of three separate validators, so the deprecated spellings
#' are translated once and every exported function inherits the same messages.
#'
#' @inheritParams validateSourceArgs
#' @param multiplexed The deprecated spelling of `multi`.
#'
#' @return A list with the normalised `manifest` and `multi`.
#'
#' @noRd
prepareSource <- function(manifest, path, multi, multiplexed = NULL,
                          call = parent.frame()) {
    multi <- resolveMultiArg(multi, multiplexed, call = call)
    validateSourceArgs(manifest, path, multi, call = call)

    if (!is.null(manifest)) {
        manifest <- migrateManifestCols(manifest, call = call)
        validateManifest(manifest, call = call)
    }

    list(manifest = manifest, multi = multi)
}

#' Validate a manifest data frame
#'
#' Column presence is checked before the row count, since a frame that is both
#' empty and malformed is more usefully reported as malformed.
#'
#' @param manifest The manifest to check.
#' @param arg Name of the caller's argument, for the error message.
#' @param call Environment to blame in the error.
#'
#' @return `invisible(TRUE)`, or aborts.
#'
#' @noRd
validateManifest <- function(manifest, arg = "manifest", call = parent.frame()) {
    required_cols <- c("Path", "multi")
    missing_cols  <- setdiff(required_cols, colnames(manifest))
    if (length(missing_cols) > 0) {
        cli::cli_abort(
            #== qty() pins the plural to the column count; without it the
            #== interpolated `arg` (length 1) would set the quantity instead
            "{.arg {arg}} is missing required {cli::qty(length(missing_cols))}column{?s}: {.field {missing_cols}}",
            call = call
        )
    }
    if (nrow(manifest) == 0) {
        cli::cli_abort("{.arg {arg}} has no rows.", call = call)
    }
    if (!is.logical(manifest$multi)) {
        cli::cli_abort(
            "{.field {paste0(arg, '$multi')}} must be logical ({.val {TRUE}}/{.val {FALSE}}).",
            call = call
        )
    }
    invisible(TRUE)
}

#' Sample name implied by a Cell Ranger `outs/` path
#'
#' The folder *above* `outs`, which is the run directory Cell Ranger was
#' pointed at and therefore the name a user recognises.
#'
#' @param outs_dir Path to an `outs/` directory, trailing slashes already
#'   stripped.
#'
#' @return A length-one character vector.
#'
#' @noRd
identFromOutsDir <- function(outs_dir) {
    basename(dirname(outs_dir))
}

#== Where `cellranger multi` puts each demuxed sample's filtered matrix,
#== relative to that sample's directory under `per_sample_outs/`. The H5 sits
#== inside a `count/` folder alongside the BAM and the analysis outputs -- it
#== has since `multi` shipped in Cell Ranger 6. The bare name is kept as a
#== fallback for pipelines that stage or flatten Cell Ranger output, mirroring
#== rawH5Candidates() below.
.multi_filtered_h5 <- c(
    file.path("count", "sample_filtered_feature_bc_matrix.h5"),
    "sample_filtered_feature_bc_matrix.h5"
)

#' Read one sample's filtered matrix into a Seurat object
#'
#' @param filename One or more paths relative to `outs_dir`, in search order;
#'   the first that exists is read. A vector rather than a single name because
#'   the same logical file lives in more than one place depending on which
#'   Cell Ranger pipeline wrote it.
#'
#' @noRd
readOneSample10x <- function(outs_dir, ident, min.cells, min.features,
                              filename = "filtered_feature_bc_matrix.h5",
                              on_missing = function(msg) cli::cli_warn("{msg}")) {
    candidates <- file.path(outs_dir, filename)
    found      <- candidates[file.exists(candidates)]

    if (length(found) == 0) {
        on_missing(cli::format_inline("File not found at {.path {candidates}}"))
        return(NULL)
    }
    mtx_path <- found[1]

    counts <- Seurat::Read10X_h5(filename = mtx_path)
    seurat <- Seurat::CreateSeuratObject(
        counts       = counts,
        project      = ident,
        min.cells    = min.cells,
        min.features = min.features,
        assay        = "RNA"
    )
    seurat$Ident <- ident

    n_cells    <- ncol(seurat)
    n_features <- nrow(seurat)
    cli::cli_alert_success(
        "{.val {ident}}: {n_cells} cell{?s}, {n_features} feature{?s}."
    )

    seurat
}

#' Read every demuxed sample in one `cellranger multi` run
#'
#' @noRd
readMultiplexedRun10x <- function(outs_dir, min.cells, min.features,
                                   on_missing = function(msg) cli::cli_warn("{msg}")) {
    per_sample_dir <- file.path(outs_dir, "per_sample_outs")
    individual_samples <- list.dirs(per_sample_dir, full.names = TRUE, recursive = FALSE)

    #== list.dirs() returns character(0) for both an absent and an empty
    #== directory, which would otherwise make the whole run vanish without a
    #== word. Almost always this means `multi = TRUE` was set on a
    #== `cellranger count` directory, so it goes through on_missing and gets
    #== the same fatal-in-path-mode treatment as an absent H5.
    if (length(individual_samples) == 0) {
        on_missing(cli::format_inline(
            "No demuxed sample directories found in {.path {per_sample_dir}}"
        ))
        return(list())
    }

    sample_seurats <- lapply(individual_samples, function(i) {
        ident <- basename(i)
        readOneSample10x(
            outs_dir     = i,
            ident        = ident,
            min.cells    = min.cells,
            min.features = min.features,
            filename     = .multi_filtered_h5,
            on_missing   = on_missing
        )
    })
    names(sample_seurats) <- basename(individual_samples)
    sample_seurats
}

# ---- raw (unfiltered) matrices -------------------------------------------
#
# The raw matrix is not the filtered matrix with a different name, and it does
# not sit where the filtered one does. `cellranger multi` writes a *single*
# run-level raw matrix under `outs/multi/count/`, shared by every demuxed
# sample in the pool, because cell calling happens once per GEM well and only
# afterwards are barcodes assigned to samples. So the fan-out that
# readMultiplexedRun10x() does over `per_sample_outs/` has no equivalent here:
# one run directory yields one raw matrix either way.

#' Where a run's raw feature-barcode H5 might live
#'
#' Returned in search order. `cellranger multi` runs get two candidates because
#' the run-level path only exists from Cell Ranger 6 onward, and because
#' pipelines that stage or flatten Cell Ranger output often drop the raw matrix
#' directly into `outs/`.
#'
#' @param outs_dir Path to an `outs/` directory.
#' @param multi Whether this is a `cellranger multi` run.
#'
#' @return A character vector of candidate paths, most likely first.
#'
#' @noRd
rawH5Candidates <- function(outs_dir, multi) {
    flat <- file.path(outs_dir, "raw_feature_bc_matrix.h5")
    if (multi) {
        c(file.path(outs_dir, "multi", "count", "raw_feature_bc_matrix.h5"), flat)
    } else {
        flat
    }
}

#' First raw H5 that actually exists for a run
#'
#' @inheritParams rawH5Candidates
#'
#' @return The path, or `NA_character_` if none of the candidates exist.
#'
#' @noRd
resolveRawH5 <- function(outs_dir, multi) {
    candidates <- rawH5Candidates(outs_dir, multi)
    found      <- candidates[file.exists(candidates)]
    if (length(found) == 0) NA_character_ else found[1]
}

#' Reduce a multi-modality Read10X_h5() result to gene expression
#'
#' A run with Antibody Capture, CRISPR or Multiplexing Capture features makes
#' [Seurat::Read10X_h5()] return a *list* of matrices rather than one matrix.
#' Barcode ranks are a gene-expression statistic, so that library is the one to
#' keep -- ranking on a hashtag library would call a knee off ~10 features.
#'
#' @param counts A matrix, or a named list of matrices.
#' @param ident Sample name, for the message.
#'
#' @return A single count matrix.
#'
#' @noRd
selectGexMatrix <- function(counts, ident) {
    if (!is.list(counts)) {
        return(counts)
    }

    modalities <- names(counts)
    if ("Gene Expression" %in% modalities) {
        cli::cli_alert_info(
            "{.val {ident}}: {length(counts)} feature types present ({.val {modalities}}); using {.val Gene Expression}."
        )
        return(counts[["Gene Expression"]])
    }

    cli::cli_warn(c(
        "{.val {ident}}: no {.val Gene Expression} feature type in the raw matrix.",
        "i" = "Using the first of {.val {modalities}} instead."
    ))
    counts[[1]]
}

#' Read one run's raw feature-barcode matrix
#'
#' Deliberately returns the count matrix rather than a Seurat object. A raw
#' matrix carries every barcode in the GEM well -- hundreds of thousands to
#' millions of columns, nearly all of them empty droplets -- and wrapping that
#' in a Seurat object costs memory to build metadata for droplets that are
#' about to be thrown away. Barcode ranking only needs the column sums.
#'
#' @param outs_dir Path to an `outs/` directory.
#' @param ident Sample name, used for messages.
#' @param multi Whether this is a `cellranger multi` run.
#' @param on_missing Callback invoked with a message when no raw H5 is found.
#'
#' @return A sparse count matrix (features x barcodes), or `NULL` if no raw H5
#'   was found.
#'
#' @noRd
readOneRaw10x <- function(outs_dir, ident, multi,
                          on_missing = function(msg) cli::cli_warn("{msg}")) {
    mtx_path <- resolveRawH5(outs_dir, multi)

    if (is.na(mtx_path)) {
        candidates <- rawH5Candidates(outs_dir, multi)
        on_missing(cli::format_inline(
            "Raw matrix for {.val {ident}} not found at {.path {candidates}}"
        ))
        return(NULL)
    }

    counts <- Seurat::Read10X_h5(filename = mtx_path)
    counts <- selectGexMatrix(counts, ident)

    n_barcodes <- ncol(counts)
    n_features <- nrow(counts)
    cli::cli_alert_success(
        "{.val {ident}}: {n_barcodes} raw barcode{?s}, {n_features} feature{?s}."
    )

    counts
}

#' Enumerate the runs a manifest or path describes, without reading anything
#'
#' Raw matrices are too large to hold several of at once, so the barcode rank
#' plotter needs the work list up front -- one entry per run -- and reads each
#' matrix only when it reaches it.
#'
#' @inheritParams validateSourceArgs
#'
#' @return A data frame with columns `ident`, `outs_dir` and `multi`, one row
#'   per run, with `ident` made unique.
#'
#' @noRd
enumerateRawRuns <- function(manifest, path, multi) {
    if (!is.null(path)) {
        outs_dirs <- sub("/+$", "", path)
        flags     <- multi
    } else {
        outs_dirs <- sub("/+$", "", as.character(manifest$Path))
        flags     <- manifest$multi
    }

    idents <- identFromOutsDir(outs_dirs)
    if (anyDuplicated(idents)) {
        dupes <- unique(idents[duplicated(idents)])
        cli::cli_warn(c(
            "Duplicate run name{?s} in the manifest: {.val {dupes}}.",
            "i" = "Names are taken from the folder above {.path outs}; suffixes were added to keep them unique."
        ))
        idents <- make.unique(idents, sep = "_")
    }

    data.frame(
        ident    = idents,
        outs_dir = outs_dirs,
        multi    = flags,
        stringsAsFactors = FALSE
    )
}
