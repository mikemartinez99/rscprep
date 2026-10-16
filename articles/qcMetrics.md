# Per-cell QC metrics and outlier flagging

``` r

library(rscprep)
```

Two functions, in order.
[`diagnoseQC()`](https://mikemartinez99.github.io/rscprep/reference/diagnoseQC.md)
calculates the per-cell metrics that filtering decisions get made from;
[`flagQC()`](https://mikemartinez99.github.io/rscprep/reference/flagQC.md)
marks the cells those metrics say are outliers. Neither one removes
anything. The motivation here is that you can look at what a threshold
would remove before you actually remove it.

This vignette picks up where
[`mergeSeurat()`](https://mikemartinez99.github.io/rscprep/articles/mergeSeurat.md)
leaves off:

``` r

seurat_list <- seuratFrom10x(manifest = manifest)
seurat_list <- includeMetadata(seurat_list, sample_sheet)
seurat      <- mergeSeurat(seurat_list)
```

Both functions want that single merged object rather than the list, so
that one set of thresholds is chosen from one view of the data.

## Calculating the metrics

``` r

seurat <- diagnoseQC(seurat, organism = "mmu")
#> ℹ Diagnosing 1200 cells in assay "RNA" using "mmu" gene patterns.
#> ℹ mito: 8 features matched "^mt".
#> ℹ ribo: 39 features matched "^Rp[sl]".
#> ✔ Added 3 QC metrics to 1200 cells.
#> • pctMito: median 0.95%, max 19.59%
#> • pctRibo: median 4.59%, max 5.91%
#> • log10GenesPerUMI: median 0.764
```

Three columns land in the metadata, alongside the `nCount_RNA` and
`nFeature_RNA` that Seurat already wrote:

``` r

head(seurat@meta.data[, c("nCount_RNA", "nFeature_RNA",
                          "pctMito", "pctRibo", "log10GenesPerUMI")])
#>                nCount_RNA nFeature_RNA pctMito pctRibo log10GenesPerUMI
#> Sample_1_Cell1        443          318  15.124   3.612            0.946
#> Sample_1_Cell2        464          326  12.500   3.879            0.943
#> Sample_1_Cell3        444          294  19.595   2.703            0.932
#> Sample_1_Cell4        454          316  15.639   3.304            0.941
```

- **`pctMito`**: share of a cell’s UMIs coming from mitochondrial genes.
  High means the cell lost cytoplasmic RNA through a ruptured membrane,
  so this is the usual first filter. What counts as high is
  tissue-dependent: 5–10% is a common starting point, but cardiac,
  muscle and kidney are legitimately mitochondria-rich and a fixed
  cutoff throws away real cells.
- **`pctRibo`**: share from ribosomal protein genes. Rarely a hard
  filter, because it varies strongly by cell type. Useful as a
  comparison between samples that ought to look alike: a sharp shift
  points at something technical.
- **`log10GenesPerUMI`**: library complexity,
  `log10(nFeature) / log10(nCount)`. Between 0 and 1 for real cells; low
  means a library dominated by few transcripts. Around 0.8 is a
  conventional floor.

The percentages come from matching gene **symbols** against
organism-specific patterns, so `organism` has to match the data. The
package currently supports mouse, rat, human, and zebrafish.

| code    | organism  | mito   | ribo      |
|---------|-----------|--------|-----------|
| `"hsa"` | human     | `^MT-` | `^RP[SL]` |
| `"mmu"` | mouse     | `^mt`  | `^Rp[sl]` |
| `"rno"` | rat       | `^Mt-` | `^Rp[sl]` |
| `"dre"` | zebrafish | `^mt-` | `^rp[sl]` |

An organism not listed can borrow the code of one whose convention it
shares.

## Cell cycle, optionally

`doCellCycle = TRUE` also hands the object to
[`runCellCycle()`](https://mikemartinez99.github.io/rscprep/reference/runCellCycle.md),
adding `S.Score`, `G2M.Score` and `Phase`:

``` r

seurat <- diagnoseQC(seurat, organism = "mmu", doCellCycle = TRUE)
```

It is off by default because, unlike the metrics above, scoring assumes
normalized data, and will run
[`Seurat::NormalizeData()`](https://satijalab.org/seurat/reference/NormalizeData.html)
on the object itself, with a warning, if it arrives raw.

## Flagging outliers

[`flagQC()`](https://mikemartinez99.github.io/rscprep/reference/flagQC.md)
marks cells that sit an unusual distance from the rest, using the median
absolute deviation rather than a fixed cutoff:

``` r

seurat <- flagQC(seurat, metrics = c(nCount_RNA = "lower", pctMito = "higher"))
```

Each entry is `column = bound`, and the bound is per metric because
which outlier tail matters depends on the metric. `nCount_RNA` wants
`"lower"`: the low tail is empty droplets and dying cells, while the
high tail is mostly real, well-captured cells. `pctMito` wants
`"higher"`, since high is bad and low is not. `"both"` is the default,
and is what you want when doublets are the concern.

A MAD-based threshold adapts to the data, which a hard cutoff cannot: a
cutoff carried between experiments discards good cells from a shallow
run and keeps bad ones from a deep one. The cost is that it assumes most
of the cells it is estimated from are good. On a badly failed sample it
will call the failure normal and flag the survivors.

## Estimating thresholds per sample

This is the argument that matters most. By default
[`flagQC()`](https://mikemartinez99.github.io/rscprep/reference/flagQC.md)
estimates one threshold across every cell in the object. If your samples
differ in sequencing depth (and they usually do), that single threshold
describes none of them.

Here are three mouse samples at median depths of roughly 2,500, 6,800
and 25,000 UMIs. Globally:

``` r

global <- flagQC(seurat, c(nCount_RNA = "lower", pctMito = "higher"))
#> ℹ Flagging outliers in 2 metrics at 3 MADs from the median (log2 scale), across all cells.
#> ℹ nCount_RNA: 0/1200 cells flagged (0%).
#> • lower: < 87.67 -- 0 cells
#> ℹ pctMito: 68/1200 cells flagged (5.67%).
#> • upper: > 1.386 -- 68 cells
```

The `nCount_RNA` threshold comes out at **87 UMIs and catches nothing**.
Pooling three depth tiers inflates the spread so far that the cutoff
falls below every real cell in the experiment, including the dying ones
it exists to find.

Now per sample:

``` r

seurat <- flagQC(
  seurat,
  metrics  = c(nCount_RNA = "lower", pctMito = "higher"),
  group_by = "Ident"
)
#> ℹ Flagging outliers in 2 metrics at 3 MADs from the median (log2 scale), within 3 groups of Ident.
#> ℹ nCount_RNA: 64/1200 cells flagged (5.33%) across 3 groups.
#> • lower: < 2400-24910 -- 64 cells
#> ℹ pctMito: 62/1200 cells flagged (5.17%) across 3 groups.
#> • upper: > 1.182-1.858 -- 62 cells
```

Three thresholds instead of one (2,400, 6,512 and 24,910), and the
low-count cells in every sample get found. The console prints the
*range* across groups rather than a line per group, so a twenty-sample
object stays readable; the per-group numbers are in the table.

`"Ident"` is the sample name
[`seuratFrom10x()`](https://mikemartinez99.github.io/rscprep/reference/seuratFrom10x.md)
sets, but any discrete metadata column works: a capture batch, a 10x
lane, a tissue. Group on the technical variable you want the threshold
to absorb.

Note that `pctMito` barely changed (68 flagged globally, 62 per sample).
It is a ratio, so it does not scale with depth. Grouping matters for
depth-dependent metrics; for ratio metrics it is close to a no-op. That
is worth checking rather than assuming.

## Reading the summary

Every call records what it did.
[`getOutlierSummary()`](https://mikemartinez99.github.io/rscprep/reference/getOutlierSummary.md)
reads it back:

``` r

qc <- getOutlierSummary(seurat)
qc[, c("metric", "group", "lower", "upper", "n_flagged", "pct_flagged")]
#>       metric    group     lower    upper n_flagged pct_flagged
#> 1 nCount_RNA Sample_1  2400.429       NA        22         5.5
#> 2 nCount_RNA Sample_2  6512.113       NA        12         3.0
#> 3 nCount_RNA Sample_3 24910.163       NA        30         7.5
#> 4    pctMito Sample_1        NA 1.857964        20         5.0
#> 5    pctMito Sample_2        NA 1.431274        12         3.0
#> 6    pctMito Sample_3        NA 1.181886        30         7.5
```

One row per metric per group. Two things to notice:

- Thresholds are on the **metric’s own scale**, even though the estimate
  was made on a log scale. They can be compared against the column
  directly.
- The unenforced tail is `NA`, not an unapplied number. `nCount_RNA` was
  flagged `"lower"` only, so its `upper` is blank rather than a cutoff
  nothing was tested against. The same applies to `n_lower` / `n_upper`.

Which sample is losing the most cells?

``` r

qc[order(-qc$pct_flagged), c("metric", "group", "n_flagged", "pct_flagged")]
#>       metric    group n_flagged pct_flagged
#> 3 nCount_RNA Sample_3        30         7.5
#> 6    pctMito Sample_3        30         7.5
#> 1 nCount_RNA Sample_1        22         5.5
#> 4    pctMito Sample_1        20         5.0
#> 2 nCount_RNA Sample_2        12         3.0
#> 5    pctMito Sample_2        12         3.0
```

Repeated calls accumulate into this one table, so an object flagged on
counts globally and then on mitochondrial content per sample yields a
single table describing both. It is a record written at flag time, not a
live view. Cell counts describe the object as it stood then, which is
the point, since being able to say afterwards how many cells were
dropped and why is the reason to keep it.

## Groups too small to estimate from

Per-group estimation creates a failure mode global estimation does not
have. A group of two cells has no usable spread, and a group where over
half the cells share one value has a MAD of zero. Both are skipped
rather than given a meaningless threshold:

``` r

seurat <- flagQC(seurat, c(nCount_RNA = "lower"), group_by = "Ident")
#> Warning:
#> ! nCount_RNA: no threshold for 1 of 3 Ident groups.
#> ℹ Skipped: tiny (fewer than 3 cells)
#> ✖ Cells in that group are left as NA, not FALSE.
#> ℹ nCount_RNA: 0/58 cells flagged (0%) across 2 groups.
```

One warning naming every skipped group, not one per group. Skipped
groups still get a row in the summary table, with `NA` thresholds and an
`NA` `column`, so nothing silently disappears.

That `NA`-not-`FALSE` distinction runs throughout. A cell is left `NA`
when it was never judged (missing metric, missing `group_by` value, or a
skipped group), so an unjudged cell is never mistaken for one that
passed. It matters at the filtering step below.

## Actually filtering

[`flagQC()`](https://mikemartinez99.github.io/rscprep/reference/flagQC.md)
writes logical columns and stops. Filtering is yours:

``` r

flags <- grep("_outliers$", colnames(seurat@meta.data), value = TRUE)
flags
#> [1] "nCount_RNA_3_lower_outliers" "pctMito_3_higher_outliers"

# `%in% TRUE` treats an unjudged NA as "not an outlier"; use `isTRUE`-style
# logic deliberately rather than letting NA propagate through `|`
keep     <- !Reduce(`|`, lapply(seurat@meta.data[flags], function(x) x %in% TRUE))
filtered <- seurat[, which(keep)]

table(seurat$Ident, keep)
#>           keep
#>            FALSE TRUE
#>   Sample_1    22  378
#>   Sample_2    12  388
#>   Sample_3    30  370
```

Cross-tabulate before you subset. A sample losing far more than its
neighbours is worth understanding rather than filtering: it may have
failed, or it may be a genuinely different tissue where the thresholds
do not transfer.

You do not have to trust every flag. The columns are independent, so
filter on the subset you believe:

``` r

filtered <- seurat[, which(!seurat$pctMito_3_higher_outliers %in% TRUE)]
```

## Tuning

``` r

# more aggressive
seurat <- flagQC(seurat, c(pctMito = "higher"), nmads = 2.5, group_by = "Ident")

# keep all but the extremes
seurat <- flagQC(seurat, c(pctMito = "higher"), nmads = 5, group_by = "Ident")

# a metric that is already symmetric does not want the log transform
seurat <- flagQC(seurat, c(log10GenesPerUMI = "lower"), log = FALSE)
```

`nmads` is in the column name, so `nmads = 3` and `nmads = 5` produce
separate columns and can be compared side by side. `group_by` is **not**
in the name, so re-running one metric with different grouping overwrites
the previous flags. The summary table is updated to match, so it never
describes a column that is no longer there.

`log = TRUE` is the default and is right for count metrics, which are
right-skewed enough that the median and MAD are misleading on the raw
scale. Turn it off for anything already roughly symmetric.

## Where this leaves you

``` r

seurat <- mergeSeurat(seuratFrom10x(manifest = manifest))
seurat <- diagnoseQC(seurat, organism = "mmu")
seurat <- flagQC(
  seurat,
  metrics  = c(nCount_RNA = "lower", nFeature_RNA = "lower", pctMito = "higher"),
  group_by = "Ident"
)

getOutlierSummary(seurat)          # what would be cut, and from where
write.csv(getOutlierSummary(seurat), "qc/outlier_summary.csv", row.names = FALSE)
```

Nothing has been removed at this point. The metrics are on the object,
the flags are on the object, and the record of how the thresholds were
chosen is on the object. Deciding what to drop is the next step, and it
is a judgement call, which is why these functions stop short of making
it.
