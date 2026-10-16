#' 1k mouse heart cells, as a Seurat object
#'
#' A cut-down copy of the 10x Genomics "1k Heart Cells from an E18 Mouse (v3
#' chemistry)" demo run, already read through [seuratFrom10x()] so examples
#' can start from a Seurat object. The same data ships as a Cell Ranger
#' `outs/` directory for the path-based readers -- see *Bundled Cell Ranger
#' output* below.
#'
#' Every cell Cell Ranger called is kept, but only 2,092 of the 31,053 genes:
#' the 2,000 detected in the most cells, plus every mitochondrial (`mt-`),
#' ribosomal protein (`Rps`/`Rpl`) and cell-cycle gene. Per-cell totals are
#' therefore about 70% of the originals, and `pctMito` / `pctRibo` from
#' [diagnoseQC()] read higher than they would on the full gene set. That is
#' fine for showing how the functions work, and wrong for drawing conclusions
#' about the tissue.
#'
#' The object was built with `min.cells = 0, min.features = 0`. With only
#' ~2,000 genes kept, the [seuratFrom10x()] defaults would drop about 50 cells
#' for an artefact of the subset rather than for anything wrong with them.
#' Reading the bundled files with the defaults gives 957 cells.
#'
#' @section Bundled Cell Ranger output:
#' `system.file("extdata", "heart_1k_v3", "outs", package = "rscprep")` is a
#' `cellranger count` `outs/` directory holding two H5 files over the same
#' genes:
#' \describe{
#'   \item{`filtered_feature_bc_matrix.h5`}{The 1,011 called cells -- what
#'     [seuratFrom10x()] reads.}
#'   \item{`raw_feature_bc_matrix.h5`}{Every barcode with at least 25 UMIs in
#'     the original run, about 20k of 6.8M -- what [rawMatrixFrom10x()] and
#'     [plotBarcodeRanks()] read. Only the lowest-ranked barcodes were
#'     dropped, so the barcodes that remain keep their true ranks.}
#' }
#'
#' @format A [SeuratObject::Seurat] object: 2,092 genes x 1,011 cells, one
#'   `RNA` assay of raw counts, with `orig.ident` and `Ident` both
#'   `"heart_1k_v3"`. Rows are mouse gene symbols.
#'
#' @source 10x Genomics, Cell Ranger 3.0.0, mm10 reference:
#'   <https://www.10xgenomics.com/datasets/1-k-heart-cells-from-an-e-18-mouse-v-3-chemistry-3-standard-3-0-0>.
#'   The full FASTQs (6.4 GB) are downloadable from the same page. The script
#'   that builds this subset is `data-raw/heart_1k_v3.R` in the package source.
#'
#' @examples
#' heart1k
#'
#' # The same data, read from the bundled Cell Ranger output
#' outs <- system.file("extdata", "heart_1k_v3", "outs", package = "rscprep")
#' seurat_list <- seuratFrom10x(path = outs, multi = FALSE)
"heart1k"
