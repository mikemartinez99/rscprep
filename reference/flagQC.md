# Flag Outlier Cells by Median Absolute Deviation

Marks cells whose QC metrics sit an unusual distance from the rest of
the data, using the median absolute deviation (MAD) rather than a fixed
cutoff. For each metric a centre and spread are estimated, a cell is
flagged if it falls more than `nmads` MADs from the centre, and the
result is written back as a logical metadata column. Nothing is subset
or removed – flagging and filtering are kept separate so the flags can
be inspected, and cross-tabulated against sample or cell type, before
any cells are dropped.

## Usage

``` r
flagQC(
  obj,
  metrics = c(nFeature_RNA = "both", nCount_RNA = "both"),
  nmads = 3,
  bound = c("both", "lower", "higher"),
  log = TRUE,
  group_by = NULL
)
```

## Arguments

- obj:

  A single
  [Seurat](https://satijalab.org/seurat/reference/reexports.html) object
  carrying the metrics to be flagged, typically the output of
  [`diagnoseQC()`](https://mikemartinez99.github.io/rscprep/reference/diagnoseQC.md).

- metrics:

  Metadata columns to flag on. Either an unnamed character vector of
  column names, which all take the function-level `bound`, or a named
  vector of `column = bound` pairs, which set the bound per metric. The
  two forms can be mixed: `c("nCount_RNA", pctMito = "higher")` flags
  counts using `bound` and mitochondrial content on the high side only.

- nmads:

  Number of MADs from the median at which a cell becomes an outlier.
  Must be a single positive number. `3` (the default) is the
  conventional choice; `2.5` is more aggressive, `5` keeps all but the
  extremes.

- bound:

  Default tail to flag, applied to every metric not given its own bound
  in `metrics`. One of `"both"` (the default), `"lower"` or `"higher"`.
  See *Choosing a bound*.

- log:

  Estimate the median and MAD on a `log2` scale. `TRUE` by default,
  which is what count metrics want. See *The log transform*.

- group_by:

  Discrete metadata column within which to estimate thresholds, so that
  each cell is judged against its own group rather than against the
  whole object. `NULL` (default) estimates one global threshold per
  metric; `"Ident"` estimates one per sample. See *Where the thresholds
  come from*.

## Value

The input object, with one logical metadata column added per metric,
named `<metric>_<nmads>_<bound>_outliers`. `TRUE` marks an outlier;
cells whose metric is `NA`, whose `group_by` value is `NA`, or whose
group was skipped stay `NA`. The column name does not record `group_by`,
so re-running the same metric with different grouping overwrites the
previous flags rather than sitting alongside them.

A summary is stored in `obj@misc$flagQC` as a data frame, one row per
metric per group, with columns:

- `metric`, `bound`, `nmads`, `scale`:

  What was asked for, and whether the estimate was made on the `"log2"`
  or `"identity"` scale.

- `group_by`, `group`:

  The column grouped on and the level this row describes, both `NA` when
  thresholds were estimated globally.

- `center`, `lower`, `upper`:

  The median and the two thresholds, on the metric's own scale. The
  threshold for a tail that was not flagged is `NA`, not an unenforced
  number.

- `n_lower`, `n_upper`:

  Cells flagged by each threshold, again `NA` for a tail that was not
  flagged. The two are disjoint, so where both apply they sum to
  `n_flagged`.

- `n_flagged`, `pct_flagged`, `n_missing`, `n_cells`:

  Totals within the group, where `n_missing` counts cells the metric was
  `NA` for and so could not be judged either way.

- `column`:

  The metadata column written, or `NA` for a group that was skipped.

Re-running on the same object appends to that table rather than
replacing it, so it stays a record of every flag column present. Rows
describing a column that is being overwritten are replaced.

## Details

MAD-based thresholds adapt to the data. A hard cutoff carried between
experiments discards good cells from a low-depth run and keeps bad ones
from a deep run; a threshold derived from the data's own distribution
does neither. The cost is that it assumes most cells it is estimated
from are good – the median and MAD are computed over all of them – so on
a badly failed sample it will happily call the failure normal and flag
the survivors.

## Where the thresholds come from

`group_by` decides what "the rest of the data" means, and it is the most
consequential argument here:

- `NULL` (default):

  One global threshold per metric, estimated across every cell in the
  object. Right for a single sample, and for a merged object whose
  samples were captured and sequenced together.

- `"Ident"`:

  A separate threshold per sample –
  [`seuratFrom10x()`](https://mikemartinez99.github.io/rscprep/reference/seuratFrom10x.md)
  sets `Ident` to the sample name. Right whenever sequencing depth or
  cell recovery varies between samples, which is the usual case. A
  global cutoff on such an object flags most of the shallowest sample
  and almost none of the deepest, regardless of whether either actually
  failed.

- any other metadata column:

  A separate threshold per level of that column: a capture batch, a 10x
  lane, a tissue. Group by the technical variable you want the threshold
  to absorb.

Grouping is not free. Each group is judged only against itself, so a
uniformly bad sample looks fine from the inside and a group with few
cells gives an unstable estimate. Groups with fewer than three cells
carrying the metric, or where over half the cells share one value, are
skipped rather than given a meaningless threshold; they are named in a
warning and appear in the summary table with `NA` thresholds. The two
views are worth comparing – if per-sample and global flags disagree
sharply, the disagreement is itself the finding.

## Choosing a bound

Which tail matters depends on the metric, so the bound is set per metric
rather than globally:

- `"lower"`:

  Only cells below the threshold. The right choice for `nFeature_RNA`
  and `nCount_RNA`, where the low tail is empty droplets and dying cells
  but the high tail is largely real, well-captured cells.

- `"higher"`:

  Only cells above it. The right choice for `pctMito`, where high is bad
  and low is not.

- `"both"`:

  Either tail. Useful for count metrics when doublets are the concern,
  since the high tail is where they sit.

## The log transform

Count metrics are right-skewed, which drags the median up and inflates
the MAD, so `log = TRUE` (the default) estimates the centre and spread
on a `log2` scale. Thresholds are transformed back before they are
applied, reported and stored, so the values in the summary table are
always on the metric's own scale and can be compared against it
directly. If a metric contains zeroes the transform becomes
`log2(x + 1)` and the back-transform compensates; negative values cannot
be log-transformed at all, and are an error rather than a silent `NaN`.

Set `log = FALSE` for metrics that are already roughly symmetric, such
as `log10GenesPerUMI` or a percentage that does not pile up against
zero.

## See also

[`diagnoseQC()`](https://mikemartinez99.github.io/rscprep/reference/diagnoseQC.md),
which calculates the metrics this flags on.

## Examples

``` r
if (FALSE) { # \dontrun{
seurat <- diagnoseQC(seurat, organism = "mmu")

# one global threshold per metric
seurat <- flagQC(seurat, c(nCount_RNA = "lower", pctMito = "higher"))

# the same, but each sample judged against itself
seurat <- flagQC(
    seurat,
    metrics  = c(nCount_RNA = "lower", nFeature_RNA = "lower", pctMito = "higher"),
    group_by = "Ident"
)

# per-sample thresholds and counts
seurat@misc$flagQC[, c("metric", "group", "lower", "upper", "n_flagged")]

# which samples lose the most cells?
qc <- seurat@misc$flagQC
with(qc, tapply(pct_flagged, group, max))

# group by a technical variable instead
seurat <- flagQC(seurat, "nCount_RNA", bound = "lower", group_by = "capture_batch")

# nothing has been removed yet -- look first
table(seurat$Ident, seurat$pctMito_3_higher_outliers)

# then filter on whichever flags you trust
flags <- grep("_outliers$", colnames(seurat@meta.data), value = TRUE)
keep  <- !Reduce(`|`, seurat@meta.data[flags])
seurat <- seurat[, which(keep)]
} # }
```
