#' Read 10x Genomics Gene Expression Sample(s) into Seurat Objects
#'
#' Reads one or more samples' filtered feature barcode matrix H5 files from
#' 10x Genomics Cell Ranger output and creates one
#' [Seurat][Seurat::CreateSeuratObject] object per sample. Supports both
#' standard `cellranger count` output (one sample per run) and multiplexed
#' `cellranger multi` output (multiple demuxed samples under
#' `outs/per_sample_outs/`), specified either as a single path or as a
#' manifest describing multiple runs.
#'
#' Exactly one of `manifest` or `path` must be supplied.
#'
#' When `path` is supplied, it must point to a Cell Ranger `outs/` directory.
#' Set `multiplexed = TRUE` if that `outs/` directory contains
#' `per_sample_outs/<sample>/sample_filtered_feature_bc_matrix.h5` for one or
#' more demuxed samples, or `multiplexed = FALSE` if it directly contains
#' `filtered_feature_bc_matrix.h5` for a single sample.
#'
#' When `manifest` is supplied, it must be a data frame with a `Path` column
#' (one `outs/` directory per row, following the same convention as `path`
#' above) and a logical `multiplexed` column indicating which rows are
#' multiplexed runs. Rows can mix multiplexed and non-multiplexed runs.
#'
#' @param manifest Data frame with columns `Path` (character; path to each
#'   run's `outs/` directory) and `multiplexed` (logical). Mutually exclusive
#'   with `path`.
#' @param path Character string. Path to a single Cell Ranger `outs/`
#'   directory. Mutually exclusive with `manifest`.
#' @param multiplexed Logical. Whether `path` refers to a multiplexed
#'   (`cellranger multi`) run. Required when `path` is supplied; ignored when
#'   `manifest` is supplied, since multiplexing status is read per-row from
#'   `manifest$multiplexed` instead.
#' @param min.cells Integer. Include features detected in at least this many
#'   cells. Passed to [Seurat::CreateSeuratObject()].
#' @param min.features Integer. Include cells where at least this many
#'   features are detected. Passed to [Seurat::CreateSeuratObject()].
#'
#' @return A named list of [Seurat][Seurat::CreateSeuratObject] objects, one
#'   per sample, each with the `RNA` assay and an `Ident` metadata column set
#'   to the sample name. List names match sample identifiers. Samples whose
#'   expected H5 file is missing are dropped from the manifest-mode result
#'   with a warning; in `path` mode, a missing file raises an error instead.
#'
#' @examples
#' \dontrun{
#' # Single non-multiplexed sample
#' seurat_list <- seuratFrom10x(
#'   path        = "/path/to/Sample_1/outs",
#'   multiplexed = FALSE
#' )
#'
#' # Single multiplexed run (multiple demuxed samples)
#' seurat_list <- seuratFrom10x(
#'   path        = "/path/to/Pool1/outs",
#'   multiplexed = TRUE
#' )
#'
#' # Manifest describing multiple runs, mixed multiplexed status
#' manifest <- data.frame(
#'   Path        = c("/path/to/Sample_1/outs", "/path/to/Pool1/outs"),
#'   multiplexed = c(FALSE, TRUE)
#' )
#' seurat_list <- seuratFrom10x(manifest = manifest, min.cells = 3, min.features = 200)
#' }
#'
#' @importFrom Seurat Read10X_h5 CreateSeuratObject
#' @export
seuratFrom10x <- function(manifest = NULL, path = NULL, multiplexed = NULL,
                           min.cells = 10, min.features = 50) {

    validateSourceArgs(manifest, path, multiplexed)

    # A missing H5 is fatal in path mode -- there is no other sample to fall
    # back on -- so the injected handler aborts instead of warning.
    fatal <- function(msg) cli::cli_abort("{msg}", call = NULL)

    if (!is.null(path)) {
        path <- sub("/+$", "", path)
        if (multiplexed) {
            seurat_list <- readMultiplexedRun10x(path, min.cells, min.features,
                                                  on_missing = fatal)
        } else {
            ident <- identFromOutsDir(path)
            seurat <- readOneSample10x(path, ident, min.cells, min.features,
                                        on_missing = fatal)
            seurat_list <- setNames(list(seurat), ident)
        }
        cli::cli_alert_success("Loaded {length(seurat_list)} sample{?s}.")
        return(seurat_list)
    }

    validateManifest(manifest)

    seurat_list <- lapply(seq_len(nrow(manifest)), function(row_i) {
        outs_dir <- sub("/+$", "", manifest$Path[row_i])
        if (manifest$multiplexed[row_i]) {
            readMultiplexedRun10x(outs_dir, min.cells, min.features)
        } else {
            ident <- identFromOutsDir(outs_dir)
            setNames(list(readOneSample10x(outs_dir, ident, min.cells, min.features)), ident)
        }
    })

    seurat_list <- unlist(seurat_list, recursive = FALSE)
    seurat_list <- seurat_list[!sapply(seurat_list, is.null)]

    n_rows <- nrow(manifest)
    cli::cli_alert_success(
        "Loaded {length(seurat_list)} sample{?s} from {n_rows} manifest row{?s}."
    )

    seurat_list
}
