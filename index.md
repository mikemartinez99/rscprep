# rscprep

> \[!WARNING\] **🚧 In development:** \*\*Vignettes and roxygen
> headers/docs are currently heavily AI written with the intent of
> getting a first-pass out. Once major sections of the package are
> finished, I will manually update these so language is clear. Testing
> is currently being done on an in-house dataset.

## Roadmap

### QC & outlier detection

**Note: some of these functions are also present in `GDSCTools`, but are
expanded in this package (better edge-case handling, better logging,
better conventions, etc…), The intent will be to remove those from
`GDSCTools` once this package is in a more mature state.**

Per-sample QC metrics calculation (nCount, nFeature, percent.mt,
percent.ribo, log10GenesPerUMI / complexity) —
[`diagnoseQC()`](https://mikemartinez99.github.io/rscprep/reference/diagnoseQC.md)

MAD-based outlier detection per metric, per sample —
[`flagQC()`](https://mikemartinez99.github.io/rscprep/reference/flagQC.md)

Sample-level thresholding strategy — **resolved, see below**

Comprehensive summary table output (per-sample QC stats, thresholds
applied, cells retained/removed) — `obj@misc$flagQC`, read back with
[`getOutlierSummary()`](https://mikemartinez99.github.io/rscprep/reference/getOutlierSummary.md)

Cell cycle scoring, optional, with ortholog conversion for non-human
organisms —
[`runCellCycle()`](https://mikemartinez99.github.io/rscprep/reference/runCellCycle.md),
or `diagnoseQC(doCellCycle = TRUE)`

**QC plotting functions**

### Next major step — QC visualization

Everything above stops at flagging. Nothing draws the distributions
those thresholds were chosen from, which is the part that makes a
threshold a judgement rather than a default:

Distribution plots per metric (violin and box, per sample and combined)

Scatter plots for the paired view (`nCount_RNA` against `nFeature_RNA`,
mito against library size), which is where doublets and ambient-heavy
cells separate from the main body

Threshold overlays: draw the cutoffs
[`flagQC()`](https://mikemartinez99.github.io/rscprep/reference/flagQC.md)
computed onto the distributions, and colour by flag status, so “what
would this remove” is answered before anything is removed

Cell cycle visualization — `S.Score` against `G2M.Score` coloured by
`Phase`, and phase composition per sample, to show whether cycling is
confounded with sample or condition

Visualize metrics on preliminary UMAPs (more stuff needed downstream to
do this though)

Groundwork already in place:
[`theme_prepro()`](https://mikemartinez99.github.io/rscprep/reference/theme_prepro.md)
and the
[`prepro_pal()`](https://mikemartinez99.github.io/rscprep/reference/prepro_pal.md)
/ `scale_*_prepro()` family, plus
[`plotBarcodeRanks()`](https://mikemartinez99.github.io/rscprep/reference/plotBarcodeRanks.md)
as a worked example of the plotting conventions to follow.

### Eventually

Bundled example dataset (`inst/extdata/`) so examples and vignettes
actually execute during `R CMD check`, rather than being `\dontrun{}` /
`eval = FALSE`

Manifest-building helper to auto-discover sample directories and
generate `Path`/`multi` columns from a top-level directory

Doublet detection integration ???

Ambient RNA / background correction utilities ???

Multi-sample merge/integration helpers built on top of
[`seuratFrom10x()`](https://mikemartinez99.github.io/rscprep/reference/seuratFrom10x.md)’s
output list

## Shipped (Preliminary)

Loading and assembly:

|  |  |
|----|----|
| [`seuratFrom10x()`](https://mikemartinez99.github.io/rscprep/reference/seuratFrom10x.md) | Read `cellranger count` or `cellranger multi` output into a named list of Seurat objects, by path or manifest |
| [`rawMatrixFrom10x()`](https://mikemartinez99.github.io/rscprep/reference/rawMatrixFrom10x.md) | The same, for raw (pre-cell-calling) matrices |
| [`includeMetadata()`](https://mikemartinez99.github.io/rscprep/reference/includeMetadata.md) | Attach a sample-level metadata sheet to each object |
| [`mergeSeurat()`](https://mikemartinez99.github.io/rscprep/reference/mergeSeurat.md) | Collapse the list into one object, namespacing barcodes by sample |

Quality control:

|  |  |
|----|----|
| [`diagnoseQC()`](https://mikemartinez99.github.io/rscprep/reference/diagnoseQC.md) | Per-cell mito %, ribo % and complexity, for human, mouse, rat or zebrafish |
| [`flagQC()`](https://mikemartinez99.github.io/rscprep/reference/flagQC.md) | MAD-based outlier flags, globally or within any metadata grouping |
| [`getOutlierSummary()`](https://mikemartinez99.github.io/rscprep/reference/getOutlierSummary.md) | The thresholds applied and cells flagged, per metric per group |
| [`runCellCycle()`](https://mikemartinez99.github.io/rscprep/reference/runCellCycle.md) | Cell cycle phase scoring, with ortholog conversion for non-human organisms |
| [`plotBarcodeRanks()`](https://mikemartinez99.github.io/rscprep/reference/plotBarcodeRanks.md) | Barcode rank (“knee”) plots from raw matrices |

Plotting helpers:
[`theme_prepro()`](https://mikemartinez99.github.io/rscprep/reference/theme_prepro.md),
[`prepro_pal()`](https://mikemartinez99.github.io/rscprep/reference/prepro_pal.md),
[`scale_fill_prepro()`](https://mikemartinez99.github.io/rscprep/reference/scale_fill_prepro.md),
[`scale_colour_prepro()`](https://mikemartinez99.github.io/rscprep/reference/scale_fill_prepro.md).

Five vignettes cover the workflow end to end, and the package carries a
testthat suite of ~590 assertions.
