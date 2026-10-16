# 1k mouse heart cells, as a Seurat object

A cut-down copy of the 10x Genomics "1k Heart Cells from an E18 Mouse
(v3 chemistry)" demo run, already read through
[`seuratFrom10x()`](https://mikemartinez99.github.io/rscprep/reference/seuratFrom10x.md)
so examples can start from a Seurat object. The same data ships as a
Cell Ranger `outs/` directory for the path-based readers – see *Bundled
Cell Ranger output* below.

## Usage

``` r
heart1k
```

## Format

A
[SeuratObject::Seurat](https://satijalab.github.io/seurat-object/reference/Seurat-class.html)
object: 2,092 genes x 1,011 cells, one `RNA` assay of raw counts, with
`orig.ident` and `Ident` both `"heart_1k_v3"`. Rows are mouse gene
symbols.

## Source

10x Genomics, Cell Ranger 3.0.0, mm10 reference:
<https://www.10xgenomics.com/datasets/1-k-heart-cells-from-an-e-18-mouse-v-3-chemistry-3-standard-3-0-0>.
The full FASTQs (6.4 GB) are downloadable from the same page. The script
that builds this subset is `data-raw/heart_1k_v3.R` in the package
source.

## Details

Every cell Cell Ranger called is kept, but only 2,092 of the 31,053
genes: the 2,000 detected in the most cells, plus every mitochondrial
(`mt-`), ribosomal protein (`Rps`/`Rpl`) and cell-cycle gene. Per-cell
totals are therefore about 70% of the originals, and `pctMito` /
`pctRibo` from
[`diagnoseQC()`](https://mikemartinez99.github.io/rscprep/reference/diagnoseQC.md)
read higher than they would on the full gene set. That is fine for
showing how the functions work, and wrong for drawing conclusions about
the tissue.

The object was built with `min.cells = 0, min.features = 0`. With only
~2,000 genes kept, the
[`seuratFrom10x()`](https://mikemartinez99.github.io/rscprep/reference/seuratFrom10x.md)
defaults would drop about 50 cells for an artefact of the subset rather
than for anything wrong with them. Reading the bundled files with the
defaults gives 957 cells.

## Bundled Cell Ranger output

`system.file("extdata", "heart_1k_v3", "outs", package = "rscprep")` is
a `cellranger count` `outs/` directory holding two H5 files over the
same genes:

- `filtered_feature_bc_matrix.h5`:

  The 1,011 called cells – what
  [`seuratFrom10x()`](https://mikemartinez99.github.io/rscprep/reference/seuratFrom10x.md)
  reads.

- `raw_feature_bc_matrix.h5`:

  Every barcode with at least 25 UMIs in the original run, about 20k of
  6.8M – what
  [`rawMatrixFrom10x()`](https://mikemartinez99.github.io/rscprep/reference/rawMatrixFrom10x.md)
  and
  [`plotBarcodeRanks()`](https://mikemartinez99.github.io/rscprep/reference/plotBarcodeRanks.md)
  read. Only the lowest-ranked barcodes were dropped, so the barcodes
  that remain keep their true ranks.

## Examples

``` r
heart1k
#> An object of class Seurat 
#> 2092 features across 1011 samples within 1 assay 
#> Active assay: RNA (2092 features, 0 variable features)
#>  2 layers present: counts, data

# The same data, read from the bundled Cell Ranger output
outs <- system.file("extdata", "heart_1k_v3", "outs", package = "rscprep")
seurat_list <- seuratFrom10x(path = outs, multi = FALSE)
#> ✔ "heart_1k_v3": 957 cells, 2089 features.
#> ✔ Loaded 1 sample.
```
