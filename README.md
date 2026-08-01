# scPrepR

> [!WARNING]
> **🚧 In development:**

## Roadmap

### v0.1.0 — QC & outlier detection (target: 2026-09-01)

- [ ] Per-sample QC metrics calculation (nCount, nFeature, percent.mt, percent.ribo, log10GenesPerUMI / complexity)
- [ ] QC plotting functions (violin/scatter plots per sample and combined across samples)
- [ ] MAD-based outlier detection per metric, per sample
- [ ] Sample-level thresholding strategy — **open design question, see below**
- [ ] Comprehensive summary table output (per-sample QC stats, thresholds applied, cells retained/removed)

**Open design question — sample thresholds:**
MAD-based outlier detection gives a statistical flag per cell per metric, but it doesn't by itself decide *where* the cutoff should live for a given sample — that still needs a rule. Worth deciding before implementation starts, not during: a fixed MAD multiplier (e.g., 3×MAD from median) applied uniformly across all samples, a per-sample-adjustable multiplier (some samples may run hotter on mito% for real biological reasons), or a hybrid — MAD-based defaults with an override mechanism so a user can hand-tune a specific sample without touching the others. The output table format depends on this decision, since it needs to show what threshold was actually used per sample, not just whether a cell passed or failed.

### Backlog (unscheduled)

- [ ] Vignette with a bundled example dataset (`inst/extdata/`) so examples in `seuratFrom10x()` docs actually execute during `R CMD check`, rather than being `eval = FALSE`
- [ ] Manifest-building helper to auto-discover sample directories and generate `Path`/`multiplexed` columns from a top-level directory
- [ ] Doublet detection integration
- [ ] Ambient RNA / background correction utilities
- [ ] Multi-sample merge/integration helpers built on top of `seuratFrom10x()`'s output list
