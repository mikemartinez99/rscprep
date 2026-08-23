# rscprep

> [!WARNING]
> **🚧 In development:**

## Roadmap

### v0.1.0 — QC & outlier detection (target: 2026-09-01)

- [x] Per-sample QC metrics calculation (nCount, nFeature, percent.mt, percent.ribo, log10GenesPerUMI / complexity) — `diagnoseQC()`
- [x] MAD-based outlier detection per metric, per sample — `flagQC()`
- [x] Sample-level thresholding strategy — **resolved, see below**
- [x] Comprehensive summary table output (per-sample QC stats, thresholds applied, cells retained/removed) — `obj@misc$flagQC`, read back with `getOutlierSummary()`
- [x] Cell cycle scoring, optional, with ortholog conversion for non-human organisms — `runCellCycle()`, or `diagnoseQC(doCellCycle = TRUE)`
- [ ] **QC plotting functions — the remaining piece, see next major step**

**Resolved — sample thresholds:**
The question was whether the MAD multiplier should be fixed across samples,
adjustable per sample, or a hybrid. It landed on the hybrid, split across three
arguments rather than one setting:

- `nmads` sets the multiplier, defaulting to 3 and applying uniformly.
- `bound` is per metric, not global, because which tail matters depends on the
  metric: `nCount_RNA` wants `"lower"`, `pctMito` wants `"higher"`.
- `group_by` decides *what the threshold is estimated from*. This turned out to
  matter more than the multiplier. With samples at different sequencing depths,
  a single global `nCount_RNA` cutoff can land below every real cell in the
  experiment and flag nothing at all, while per-sample estimation finds the bad
  cells in each. `group_by = "Ident"` is per sample; any discrete metadata
  column works, so a capture batch or lane can absorb the variation instead.

A sample that genuinely runs hotter on mito% is therefore handled by grouping,
not by hand-tuning a multiplier for it. The summary table records the
thresholds actually applied per metric per group, so what was cut and why stays
inspectable after the fact.

### Next major step — QC visualization

Everything above stops at flagging. Nothing draws the distributions those
thresholds were chosen from, which is the part that makes a threshold a
judgement rather than a default:

- [ ] Distribution plots per metric — violin and box, per sample and combined
- [ ] Scatter plots for the paired view (`nCount_RNA` against `nFeature_RNA`,
      mito against library size), which is where doublets and ambient-heavy
      cells separate from the main body
- [ ] Threshold overlays: draw the cutoffs `flagQC()` computed onto the
      distributions, and colour by flag status, so "what would this remove"
      is answered before anything is removed
- [ ] Cell cycle visualization — `S.Score` against `G2M.Score` coloured by
      `Phase`, and phase composition per sample, to show whether cycling is
      confounded with sample or condition

Groundwork already in place: `theme_prepro()` and the `prepro_pal()` /
`scale_*_prepro()` family, plus `plotBarcodeRanks()` as a worked example of the
plotting conventions to follow.

### Backlog (unscheduled)

- [ ] Bundled example dataset (`inst/extdata/`) so examples and vignettes
      actually execute during `R CMD check`, rather than being `\dontrun{}` /
      `eval = FALSE`
- [ ] Manifest-building helper to auto-discover sample directories and generate `Path`/`multi` columns from a top-level directory
- [ ] Doublet detection integration
- [ ] Ambient RNA / background correction utilities
- [ ] Multi-sample merge/integration helpers built on top of `seuratFrom10x()`'s output list
- [ ] Fill in the `License` field — `DESCRIPTION` still holds the `usethis`
      placeholder, which fails `R CMD check`

## Shipped

Loading and assembly:

| | |
|---|---|
| `seuratFrom10x()` | Read `cellranger count` or `cellranger multi` output into a named list of Seurat objects, by path or manifest |
| `rawMatrixFrom10x()` | The same, for raw (pre-cell-calling) matrices |
| `includeMetadata()` | Attach a sample-level metadata sheet to each object |
| `mergeSeurat()` | Collapse the list into one object, namespacing barcodes by sample |

Quality control:

| | |
|---|---|
| `diagnoseQC()` | Per-cell mito %, ribo % and complexity, for human, mouse, rat or zebrafish |
| `flagQC()` | MAD-based outlier flags, globally or within any metadata grouping |
| `getOutlierSummary()` | The thresholds applied and cells flagged, per metric per group |
| `runCellCycle()` | Cell cycle phase scoring, with ortholog conversion for non-human organisms |
| `plotBarcodeRanks()` | Barcode rank ("knee") plots from raw matrices |

Plotting helpers: `theme_prepro()`, `prepro_pal()`, `scale_fill_prepro()`,
`scale_colour_prepro()`.

Five vignettes cover the workflow end to end, and the package carries a
testthat suite of ~590 assertions.
