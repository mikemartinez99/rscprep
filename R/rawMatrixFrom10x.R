#' Read 10x Genomics Raw (Unfiltered) Feature Barcode Matrices
#'
#' The raw-matrix counterpart of [seuratFrom10x()]. Reads the *unfiltered*
#' feature barcode matrix H5 -- every barcode Cell Ranger saw, not just the
#' ones it called as cells -- for one or more runs, taking the same
#' `manifest` / `path` / `multiplexed` arguments.
#'
#' This is the input a barcode rank plot needs. Cell calling is exactly the
#' step that separates the raw matrix from the filtered one, so a plot meant to
#' show *where that cutoff landed* cannot be drawn from the filtered matrix:
#' the empty droplets that form the left-hand cliff have already been removed.
#' Use [plotBarcodeRanks()] if drawing that plot is all you want; it calls this
#' reader one run at a time and never holds more than one raw matrix in memory.
#'
#' @section What is returned, and why it is not Seurat objects:
#'
#' Sparse count matrices, not [Seurat][Seurat::CreateSeuratObject] objects. A
#' raw matrix has one column per barcode in the GEM well -- of order 10^5 to
#' 10^6, the overwhelming majority of them empty droplets -- and building a
#' Seurat object means computing and storing per-cell metadata for all of them.
#' Nothing downstream of a raw matrix wants that; barcode ranking needs only
#' the column sums. Filtered matrices, which are already cell-called, are what
#' [seuratFrom10x()] wraps in Seurat objects.
#'
#' @section Multiplexed runs yield one matrix per run, not per sample:
#'
#' This is the one place where the return shape differs from
#' [seuratFrom10x()]. `cellranger multi` writes a single run-level raw matrix
#' at `outs/multi/count/raw_feature_bc_matrix.h5`, shared by every demuxed
#' sample in the pool -- cell calling happens once for the whole GEM well, and
#' only afterwards are the called barcodes assigned to samples. There is
#' therefore no per-sample raw matrix to read, and a multiplexed row expands to
#' one entry named after the run directory rather than to N entries named after
#' the demuxed samples.
#'
#' Read that way, a barcode rank curve from a multiplexed run describes the
#' pool, which is the level at which the knee is a meaningful diagnostic.
#'
#' @inheritParams seuratFrom10x
#'
#' @return A named list of sparse count matrices (features x barcodes), one per
#'   run, named after the directory above each `outs/`. Runs whose raw H5 is
#'   missing are dropped from the manifest-mode result with a warning; in
#'   `path` mode a missing file raises an error instead.
#'
#' @examples
#' \dontrun{
#' # One run
#' raw <- rawMatrixFrom10x(path = "/path/to/Sample_1/outs", multiplexed = FALSE)
#' dim(raw[["Sample_1"]])
#'
#' # A multiplexed pool -- one matrix for the whole run, not one per sample
#' raw <- rawMatrixFrom10x(path = "/path/to/Pool1/outs", multiplexed = TRUE)
#' names(raw)
#'
#' # Several runs from a manifest. Watch memory: these are large.
#' manifest <- data.frame(
#'   Path        = c("/path/to/Sample_1/outs", "/path/to/Pool1/outs"),
#'   multiplexed = c(FALSE, TRUE)
#' )
#' raw <- rawMatrixFrom10x(manifest = manifest)
#'
#' # Column sums are the barcode totals a rank plot is built from
#' totals <- Matrix::colSums(raw[[1]])
#' }
#'
#' @seealso [plotBarcodeRanks()], which consumes these matrices one at a time;
#'   [seuratFrom10x()] for the filtered, cell-called matrices.
#'
#' @importFrom Seurat Read10X_h5
#' @export
rawMatrixFrom10x <- function(manifest = NULL, path = NULL, multiplexed = NULL) {

    validateSourceArgs(manifest, path, multiplexed)
    if (!is.null(manifest)) validateManifest(manifest)

    # A missing H5 is fatal in path mode -- there is no other run to fall back
    # on -- so the injected handler aborts instead of warning.
    fatal      <- function(msg) cli::cli_abort("{msg}", call = NULL)
    on_missing <- if (is.null(path)) function(msg) cli::cli_warn("{msg}") else fatal

    runs <- enumerateRawRuns(manifest, path, multiplexed)

    if (nrow(runs) > 1) {
        cli::cli_alert_info(
            "Reading {nrow(runs)} raw matri{?x/ces} into memory at once; each holds every barcode in its GEM well."
        )
    }

    raw_list <- lapply(seq_len(nrow(runs)), function(i) {
        readOneRaw10x(
            outs_dir    = runs$outs_dir[i],
            ident       = runs$ident[i],
            multiplexed = runs$multiplexed[i],
            on_missing  = on_missing
        )
    })
    names(raw_list) <- runs$ident

    raw_list <- raw_list[!vapply(raw_list, is.null, logical(1))]

    cli::cli_alert_success(
        "Loaded {length(raw_list)} raw matri{?x/ces} from {nrow(runs)} run{?s}."
    )

    raw_list
}
