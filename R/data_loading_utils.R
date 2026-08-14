#== Internal workers shared by seuratFrom10x(), rawMatrixFrom10x() and
#== plotBarcodeRanks().
#==
#== Two things live here. First, the argument validation both entry points
#== share -- they take the same `manifest`/`path`/`multiplexed` triple, so the
#== checks and their messages are defined once. Second, the actual file I/O.
#==
#== Every reader takes an `on_missing` callback rather than branching on a mode
#== flag: the default warns and returns NULL (manifest mode, where one absent
#== sample should not kill a batch), and `path` mode injects a handler that
#== aborts instead.

#' Validate the `manifest` / `path` / `multiplexed` triple
#'
#' @param manifest,path,multiplexed The caller's arguments, unmodified.
#' @param call Environment to blame in the error, so the message points at the
#'   exported function rather than at this helper.
#'
#' @return `invisible(TRUE)`, or aborts.
#'
#' @noRd
validateSourceArgs <- function(manifest, path, multiplexed, call = parent.frame()) {
    if (is.null(manifest) && is.null(path)) {
        cli::cli_abort("Provide either {.arg manifest} or {.arg path}.", call = call)
    }
    if (!is.null(manifest) && !is.null(path)) {
        cli::cli_abort("Provide only one of {.arg manifest} or {.arg path}, not both.",
                       call = call)
    }
    if (!is.null(path) && is.null(multiplexed)) {
        cli::cli_abort(
            "When using {.arg path}, you must also specify {.arg multiplexed} ({.val {TRUE}}/{.val {FALSE}}).",
            call = call
        )
    }
    invisible(TRUE)
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
    required_cols <- c("Path", "multiplexed")
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
    if (!is.logical(manifest$multiplexed)) {
        cli::cli_abort(
            "{.field {paste0(arg, '$multiplexed')}} must be logical ({.val {TRUE}}/{.val {FALSE}}).",
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

#' @noRd
readOneSample10x <- function(outs_dir, ident, min.cells, min.features,
                              filename = "filtered_feature_bc_matrix.h5",
                              on_missing = function(msg) cli::cli_warn("{msg}")) {
    mtx_path <- file.path(outs_dir, filename)
    if (!file.exists(mtx_path)) {
        on_missing(cli::format_inline("File not found at {.path {mtx_path}}"))
        return(NULL)
    }
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

#' @noRd
readMultiplexedRun10x <- function(outs_dir, min.cells, min.features,
                                   on_missing = function(msg) cli::cli_warn("{msg}")) {
    per_sample_dir <- file.path(outs_dir, "per_sample_outs")
    individual_samples <- list.dirs(per_sample_dir, full.names = TRUE, recursive = FALSE)

    sample_seurats <- lapply(individual_samples, function(i) {
        ident <- basename(i)
        readOneSample10x(
            outs_dir     = i,
            ident        = ident,
            min.cells    = min.cells,
            min.features = min.features,
            filename     = "sample_filtered_feature_bc_matrix.h5",
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
#' Returned in search order. Multiplexed runs get two candidates because the
#' run-level path only exists from Cell Ranger 6 onward, and because pipelines
#' that stage or flatten Cell Ranger output often drop the raw matrix directly
#' into `outs/`.
#'
#' @param outs_dir Path to an `outs/` directory.
#' @param multiplexed Whether this is a `cellranger multi` run.
#'
#' @return A character vector of candidate paths, most likely first.
#'
#' @noRd
rawH5Candidates <- function(outs_dir, multiplexed) {
    flat <- file.path(outs_dir, "raw_feature_bc_matrix.h5")
    if (multiplexed) {
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
resolveRawH5 <- function(outs_dir, multiplexed) {
    candidates <- rawH5Candidates(outs_dir, multiplexed)
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
#' @param multiplexed Whether this is a `cellranger multi` run.
#' @param on_missing Callback invoked with a message when no raw H5 is found.
#'
#' @return A sparse count matrix (features x barcodes), or `NULL` if no raw H5
#'   was found.
#'
#' @noRd
readOneRaw10x <- function(outs_dir, ident, multiplexed,
                          on_missing = function(msg) cli::cli_warn("{msg}")) {
    mtx_path <- resolveRawH5(outs_dir, multiplexed)

    if (is.na(mtx_path)) {
        candidates <- rawH5Candidates(outs_dir, multiplexed)
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
#' @return A data frame with columns `ident`, `outs_dir` and `multiplexed`, one
#'   row per run, with `ident` made unique.
#'
#' @noRd
enumerateRawRuns <- function(manifest, path, multiplexed) {
    if (!is.null(path)) {
        outs_dirs <- sub("/+$", "", path)
        flags     <- multiplexed
    } else {
        outs_dirs <- sub("/+$", "", as.character(manifest$Path))
        flags     <- manifest$multiplexed
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
        ident       = idents,
        outs_dir    = outs_dirs,
        multiplexed = flags,
        stringsAsFactors = FALSE
    )
}
