# Builds the bundled example data from the 10x Genomics "1k Heart Cells from
# an E18 mouse (v3 chemistry)" demo run (Cell Ranger 3.0.0, mm10):
# https://www.10xgenomics.com/datasets/1-k-heart-cells-from-an-e-18-mouse-v-3-chemistry-3-standard-3-0-0
#
# Writes two things:
#   * inst/extdata/heart_1k_v3/outs/{filtered,raw}_feature_bc_matrix.h5 -- a
#     cut-down `cellranger count` outs/ directory, so the path-based readers
#     (seuratFrom10x(), rawMatrixFrom10x(), plotBarcodeRanks()) have a real
#     layout to point at.
#   * data/heart1k.rda -- the filtered matrix already loaded through
#     seuratFrom10x(), for examples that start from a Seurat object.
#
# The originals are 5 MB (filtered) and 133 MB (raw), so both are cut down:
#
#   Genes    -- the 2,000 detected in the most cells, plus every mitochondrial,
#               ribosomal protein and cell-cycle gene so diagnoseQC() and
#               runCellCycle() still find what they look for. The same gene set
#               is used in both matrices, so a cell's counts are identical in
#               each and the raw matrix stays a superset of the filtered one.
#   Barcodes -- the raw matrix keeps every barcode with >= 25 UMIs (before the
#               gene subset), about 20k of the original 6.8M. Only the lowest
#               ranks are dropped, so every kept barcode keeps its true rank
#               and the knee between cells (>= 506 UMIs) and empty droplets is
#               intact. Sampling the tail instead would shift every rank.
#
# Run from the package root:  source("data-raw/heart_1k_v3.R")

library(Matrix)

base_url <- "https://cf.10xgenomics.com/samples/cell-exp/3.0.0/heart_1k_v3/"
n_top    <- 2000
min_umi  <- 25

dl_dir <- file.path(tempdir(), "heart_1k_v3")
dir.create(dl_dir, showWarnings = FALSE)
fetch <- function(kind) {
    dest <- file.path(dl_dir, sprintf("heart_1k_v3_%s_feature_bc_matrix.h5", kind))
    if (!file.exists(dest)) {
        options(timeout = max(600, getOption("timeout")))
        utils::download.file(paste0(base_url, basename(dest)), dest, mode = "wb")
    }
    dest
}

# Read10X_h5() keys rows by symbol or by Ensembl ID, never both; read twice to
# keep both columns of features.tsv. Row order is the same either way.
read_both <- function(h5) {
    m   <- Seurat::Read10X_h5(h5, use.names = TRUE)
    ids <- rownames(Seurat::Read10X_h5(h5, use.names = FALSE))
    list(counts = m, ids = stats::setNames(ids, rownames(m)))
}

filt <- read_both(fetch("filtered"))
raw  <- read_both(fetch("raw"))
stopifnot(identical(rownames(filt$counts), rownames(raw$counts)))

# ---- gene subset -----------------------------------------------------------

genes     <- rownames(filt$counts)
cc_genes  <- unlist(Seurat::cc.genes.updated.2019, use.names = FALSE)
must_keep <- genes[grepl("^mt-|^Rp[sl]", genes) | toupper(genes) %in% cc_genes]
detected  <- Matrix::rowSums(filt$counts > 0)
top       <- names(sort(detected, decreasing = TRUE))[seq_len(n_top)]
keep      <- genes[genes %in% union(top, must_keep)]   # original order

# ---- barcode subset --------------------------------------------------------

totals   <- Matrix::colSums(raw$counts)
keep_bc  <- totals >= min_umi
stopifnot(all(colnames(filt$counts) %in% colnames(raw$counts)[keep_bc]))

# ---- write -----------------------------------------------------------------

outs <- file.path("inst", "extdata", "heart_1k_v3", "outs")
dir.create(outs, recursive = TRUE, showWarnings = FALSE)

write_h5 <- function(counts, path) {
    if (file.exists(path)) file.remove(path)
    DropletUtils::write10xCounts(
        path        = path,
        x           = counts,
        gene.id     = unname(filt$ids[rownames(counts)]),
        gene.symbol = rownames(counts),
        gene.type   = "Gene Expression",
        genome      = "mm10",
        type        = "HDF5",
        version     = "3"
    )
}

write_h5(filt$counts[keep, ], file.path(outs, "filtered_feature_bc_matrix.h5"))
write_h5(raw$counts[keep, keep_bc], file.path(outs, "raw_feature_bc_matrix.h5"))

# ---- Seurat object ---------------------------------------------------------

# Built from the bundled H5 rather than the original, so the object and the
# files describe exactly the same data. The read-time filters are off: with
# only ~2,000 genes kept, the default min.features = 50 would drop ~50 cells
# for an artefact of the subset rather than for anything wrong with them, and
# the low-quality cells are what the QC examples need to show flagQC() finding.
pkgload::load_all(".")
heart1k <- seuratFrom10x(path = outs, multi = FALSE,
                         min.cells = 0, min.features = 0)[["heart_1k_v3"]]
usethis::use_data(heart1k, overwrite = TRUE, compress = "xz")
