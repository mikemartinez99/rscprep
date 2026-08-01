# Test fixtures for scPrepR.
#
# Strategy: the package's real logic is *path resolution* -- deciding which H5
# file to look for, what to name the resulting sample, and whether a missing
# file is fatal. None of that requires a real H5 file, so fixtures create the
# Cell Ranger directory layouts with empty placeholder files (the code only
# ever calls file.exists() on them) and mock Seurat::Read10X_h5() to return an
# in-memory counts matrix. Seurat::CreateSeuratObject() is NOT mocked -- it
# runs for real, so the tests also verify we hand it valid arguments.

# ---- counts matrix -------------------------------------------------------

#' Build a small sparse counts matrix.
#'
#' With `extras = TRUE` the matrix contains two deliberately marginal entries
#' used to test min.cells / min.features pass-through:
#'   * `RareGene`   -- detected in only 2 cells (dropped by min.cells >= 3)
#'   * `SparseCell` -- only 5 features detected (dropped by min.features >= 6)
#' Everything else is dense enough to survive the package defaults
#' (min.cells = 10, min.features = 50), so the default-argument case yields a
#' predictable 60 features x 20 cells.
make_counts <- function(n_genes = 60, n_cells = 20, extras = TRUE) {
    m <- matrix(
        1,
        nrow     = n_genes,
        ncol     = n_cells,
        dimnames = list(paste0("Gene", seq_len(n_genes)),
                        paste0("Cell", seq_len(n_cells)))
    )

    if (extras) {
        rare <- matrix(0, nrow = 1, ncol = ncol(m),
                       dimnames = list("RareGene", colnames(m)))
        rare[1, 1:2] <- 1
        m <- rbind(m, rare)

        sparse <- matrix(0, nrow = nrow(m), ncol = 1,
                         dimnames = list(rownames(m), "SparseCell"))
        sparse[1:5, 1] <- 1
        m <- cbind(m, sparse)
    }

    Matrix::Matrix(m, sparse = TRUE)
}

# ---- directory layouts ---------------------------------------------------

#' Create a `cellranger count` layout: <root>/<sample>/outs/filtered_...h5
#' Returns the `outs` directory path.
make_count_run <- function(root, sample = "Sample_1", h5 = TRUE) {
    outs <- file.path(root, sample, "outs")
    dir.create(outs, recursive = TRUE, showWarnings = FALSE)
    if (h5) file.create(file.path(outs, "filtered_feature_bc_matrix.h5"))
    outs
}

#' Create a `cellranger multi` layout:
#'   <root>/<run>/outs/per_sample_outs/<sample>/sample_filtered_...h5
#' `h5` may be a single logical or one logical per sample, so a run can be
#' built with some samples' H5 files deliberately absent.
#' Returns the `outs` directory path.
make_multi_run <- function(root, run = "Pool1", samples = c("Sample_A", "Sample_B"),
                           h5 = TRUE, make_per_sample_outs = TRUE) {
    outs <- file.path(root, run, "outs")
    dir.create(outs, recursive = TRUE, showWarnings = FALSE)

    if (make_per_sample_outs) {
        h5 <- rep_len(h5, length(samples))
        for (i in seq_along(samples)) {
            sdir <- file.path(outs, "per_sample_outs", samples[i])
            dir.create(sdir, recursive = TRUE, showWarnings = FALSE)
            if (h5[i]) {
                file.create(file.path(sdir, "sample_filtered_feature_bc_matrix.h5"))
            }
        }
    }

    outs
}

# ---- mocking -------------------------------------------------------------

#' Mock Seurat::Read10X_h5() for the duration of the calling test.
#'
#' Returns an environment whose `$files` field accumulates every `filename`
#' the package asked for, in call order -- this is how tests assert that the
#' correct H5 path was resolved, not merely that *a* file was read.
local_mock_read10x <- function(counts = make_counts(), .env = parent.frame()) {
    calls <- new.env(parent = emptyenv())
    calls$files <- character(0)

    testthat::local_mocked_bindings(
        Read10X_h5 = function(filename, ...) {
            calls$files <- c(calls$files, filename)
            counts
        },
        .package = "Seurat",
        .env     = .env
    )

    calls
}

#' A temp directory scoped to the calling test.
local_root <- function(.env = parent.frame()) {
    withr::local_tempdir(.local_envir = .env)
}
