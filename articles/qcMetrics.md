# Per-cell QC metrics and outlier flagging

``` r

library(rscprep)
```

Two functions, in order.
[`diagnoseQC()`](https://mikemartinez99.github.io/rscprep/reference/diagnoseQC.md)
calculates the per-cell metrics that filtering decisions get made from;
[`flagQC()`](https://mikemartinez99.github.io/rscprep/reference/flagQC.md)
marks the cells those metrics say are outliers. Neither one removes
anything – the split is deliberate, so that you can look at what a
threshold would cut before you cut it.

## The example data

Every example in this vignette runs on data bundled with rscprep: a
cut-down copy of the 10x Genomics demo run [1k Heart Cells from an E18
Mouse (v3
chemistry)](https://www.10xgenomics.com/datasets/1-k-heart-cells-from-an-e-18-mouse-v-3-chemistry-3-standard-3-0-0).
It is one sample from one `cellranger count` run (Cell Ranger 3.0.0,
mm10 reference), with all 1,011 cells Cell Ranger called.

It comes in two forms, holding the same data:

- **`heart1k`** – a Seurat object, for examples that start after
  loading.
- **A Cell Ranger `outs/` directory** – for the functions that take a
  path. It holds `filtered_feature_bc_matrix.h5` (the called cells) and
  `raw_feature_bc_matrix.h5` (every barcode that had at least 25 UMIs in
  the original run, empty droplets included):

``` r

outs <- system.file("extdata", "heart_1k_v3", "outs", package = "rscprep")
```

To keep the package small, only 2,092 of the 31,053 genes are kept: the
2,000 detected in the most cells, plus every mitochondrial, ribosomal
protein and cell-cycle gene. Each cell therefore keeps about 70% of its
UMIs, and percentages such as `pctMito` read higher than they would on
the full gene set. That is fine for showing how the functions behave and
wrong for drawing conclusions about heart tissue.
[`?heart1k`](https://mikemartinez99.github.io/rscprep/reference/heart1k.md)
has the details, and `data-raw/heart_1k_v3.R` in the package source
rebuilds it from the 10x download.

### Making several samples from one

`heart1k` is a single sample, and this vignette needs several. So the
cells are dealt at random into three pseudo-samples, and two of them are
thinned to look like shallower sequencing runs:

``` r

counts <- heart1k[["RNA"]]@counts

set.seed(1)
sample_of <- sample(c("Sample_A", "Sample_B", "Sample_C"), ncol(counts),
                    replace = TRUE)

# Fraction of each UMI count kept. Keeping each UMI with probability p
# ("binomial thinning") is what sequencing a library p times as deeply does.
depth <- c(Sample_A = 0.15, Sample_B = 0.4, Sample_C = 1)

seurat_list <- lapply(names(depth), function(s) {
  m    <- counts[, sample_of == s]
  m@x  <- as.numeric(stats::rbinom(length(m@x), m@x, depth[[s]]))
  m    <- Matrix::drop0(m)
  obj  <- Seurat::CreateSeuratObject(m, project = s,
                                     min.cells = 0, min.features = 0)
  obj$Ident <- s     # the sample-name column seuratFrom10x() would set
  obj
})
names(seurat_list) <- names(depth)

seurat_list
#> $Sample_A
#> An object of class Seurat 
#> 2092 features across 356 samples within 1 assay 
#> Active assay: RNA (2092 features, 0 variable features)
#>  1 layer present: counts
#> 
#> $Sample_B
#> An object of class Seurat 
#> 2092 features across 341 samples within 1 assay 
#> Active assay: RNA (2092 features, 0 variable features)
#>  1 layer present: counts
#> 
#> $Sample_C
#> An object of class Seurat 
#> 2092 features across 314 samples within 1 assay 
#> Active assay: RNA (2092 features, 0 variable features)
#>  1 layer present: counts
```

That list has the same shape
[`seuratFrom10x()`](https://mikemartinez99.github.io/rscprep/reference/seuratFrom10x.md)
returns for a multi-sample run: named by sample, each object carrying
its name in `Ident`. The cells are real, but the three “samples” are not
separate biology. They differ only in simulated depth.

The three depths matter here: they are what the per-sample section below
is about.

## Where this starts

This vignette picks up where
[`mergeSeurat()`](https://mikemartinez99.github.io/rscprep/articles/mergeSeurat.md)
leaves off. With your own data that is:

``` r

seurat_list <- seuratFrom10x(manifest = manifest)
seurat_list <- includeMetadata(seurat_list, sample_sheet)
seurat      <- mergeSeurat(seurat_list)
```

Here it is the pseudo-sample list:

``` r

seurat <- mergeSeurat(seurat_list)
#> ℹ Merging 3 seurat objects...
#> ✔ Successfully merged 3 seurat objects.
```

``` r

tapply(seurat$nCount_RNA, seurat$Ident, median)
#> Sample_A Sample_B Sample_C 
#>      699     2343     5348
```

Both functions want that single merged object rather than the list, so
that one set of thresholds is chosen from one view of the data.

## Calculating the metrics

``` r

seurat <- diagnoseQC(seurat, organism = "mmu")
#> ℹ Diagnosing 1011 cells in assay "RNA" using "mmu" gene patterns.
#> ℹ mito: 13 features matched "^mt".
#> ℹ ribo: 105 features matched "^Rp[sl]".
#> ✔ Added 3 QC metrics to 1011 cells.
#> • pctMito: median 9.95%, max 99.48%
#> • pctRibo: median 20.39%, max 56.36%
#> • log10GenesPerUMI: median 0.829
```

Three columns land in the metadata, alongside the `nCount_RNA` and
`nFeature_RNA` that Seurat already wrote:

``` r

head(seurat@meta.data[, c("nCount_RNA", "nFeature_RNA",
                          "pctMito", "pctRibo", "log10GenesPerUMI")])
#>                             nCount_RNA nFeature_RNA   pctMito   pctRibo
#> Sample_A_AAACCCACACCAGTAT-1       1718          647  7.101281 31.490105
#> Sample_A_AAACCCAGTGGAACAC-1       1972          782  9.178499 23.326572
#> Sample_A_AAACGAAGTCAGGTGA-1        495           25  0.000000  4.040404
#> Sample_A_AACAACCAGGTTTGAA-1       1228          525 11.889251 35.016287
#> Sample_A_AACAAGATCGGATACT-1        916          428 36.899563  8.624454
#> Sample_A_AACAAGATCGGTAAGG-1       1654          684  8.464329 25.272068
#>                             log10GenesPerUMI
#> Sample_A_AAACCCACACCAGTAT-1        0.8688977
#> Sample_A_AAACCCAGTGGAACAC-1        0.8780845
#> Sample_A_AAACGAAGTCAGGTGA-1        0.5187921
#> Sample_A_AACAACCAGGTTTGAA-1        0.8805389
#> Sample_A_AACAAGATCGGATACT-1        0.8884324
#> Sample_A_AACAAGATCGGTAAGG-1        0.8808528
```

- **`pctMito`** – share of a cell’s UMIs coming from mitochondrial
  genes. High means the cell lost cytoplasmic RNA through a ruptured
  membrane, so this is the usual first filter. What counts as high is
  tissue-dependent: 5–10% is a common starting point, but cardiac,
  muscle and kidney are legitimately mitochondria-rich and a fixed
  cutoff throws away real cells. This is heart: 50% of its cells are
  above 10%, and 30% above 20%. (The gene subset inflates these somewhat
  – see above.)
- **`pctRibo`** – share from ribosomal protein genes. Rarely a hard
  filter, because it varies strongly by cell type. Useful as a
  comparison between samples that ought to look alike: a sharp shift
  points at something technical.
- **`log10GenesPerUMI`** – library complexity,
  `log10(nFeature) / log10(nCount)`. Between 0 and 1 for real cells; low
  means a library dominated by few transcripts. Around 0.8 is a
  conventional floor.

The percentages come from matching gene **symbols** against
organism-specific patterns, so `organism` has to match the data.

## When the patterns match nothing

Get the organism wrong and both columns come back as zeroes. That is
loud rather than silent:

``` r

wrong <- diagnoseQC(seurat, organism = "hsa")   # mouse data
#> ℹ Diagnosing 1011 cells in assay "RNA" using "hsa" gene patterns.
#> Warning: ! mito pattern "^MT-" matched 0 features for organism "hsa".
#> ℹ Check that rownames are gene symbols and not Ensembl IDs --
#>   `head(rownames(seurat))` will show which.
#> ℹ pctMito will be 0 for every cell.
#> Warning: ! ribo pattern "^RP[SL]" matched 0 features for organism "hsa".
#> ℹ Check that rownames are gene symbols and not Ensembl IDs --
#>   `head(rownames(seurat))` will show which.
#> ℹ pctRibo will be 0 for every cell.
#> ✔ Added 3 QC metrics to 1011 cells.
#> • pctMito: median 0%, max 0%
#> • pctRibo: median 0%, max 0%
#> • log10GenesPerUMI: median 0.829
```

The same warning fires when the organism is right but the matrix is
still keyed by Ensembl IDs, which is the other common cause –
`ENSMUSG00000064341` matches no symbol pattern. Convert the rownames to
symbols before diagnosing.

Four organisms are supported, differing only in the capitalisation
convention their symbols follow:

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
adding `S.Score`, `G2M.Score` and `Phase`. It is not run here, because
for mouse data it converts the human marker lists through `orthogene`,
which needs network access:

``` r

seurat <- diagnoseQC(seurat, organism = "mmu", doCellCycle = TRUE)
```

It is off by default because, unlike the metrics above, scoring assumes
normalised data – and will run
[`Seurat::NormalizeData()`](https://satijalab.org/seurat/reference/NormalizeData.html)
on the object itself, with a warning, if it arrives raw. See
[`?runCellCycle`](https://mikemartinez99.github.io/rscprep/reference/runCellCycle.md)
for the detail, including the `method` argument for machines without
network access.

## Flagging outliers

[`flagQC()`](https://mikemartinez99.github.io/rscprep/reference/flagQC.md)
marks cells that sit an unusual distance from the rest, using the median
absolute deviation rather than a fixed cutoff:

``` r

flagged <- flagQC(seurat, metrics = c(nFeature_RNA = "lower", pctMito = "higher"))
#> ℹ Flagging outliers in 2 metrics at 3 MADs from the median (log2 scale), across all cells.
#> ℹ nFeature_RNA: 14/1011 cells flagged (1.38%).
#> • lower: < 10.91 -- 14 cells
#> ℹ pctMito: 0/1011 cells flagged (0%).
#> • upper: > 167.5 -- 0 cells
#> ✔ Added 2 outlier columns; summary in `seurat@misc$flagQC`.
```

Each entry is `column = bound`, and the bound is per metric because
which tail matters depends on the metric. `nFeature_RNA` wants `"lower"`
– the low tail is empty droplets and dying cells, while the high tail is
mostly real, well-captured cells. `pctMito` wants `"higher"`, since high
is bad and low is not. `"both"` is the default, and is what you want
when doublets are the concern.

A MAD-based threshold adapts to the data, which a hard cutoff cannot: a
cutoff carried between experiments discards good cells from a shallow
run and keeps bad ones from a deep one. The cost is that it assumes most
of the cells it is estimated from are good. On a badly failed sample it
will call the failure normal and flag the survivors.

### When the log scale is wrong

Look at the `pctMito` line above: an upper cutoff of 167%, which no cell
can exceed. The MAD is estimated on a log2 scale by default. That suits
count metrics, which are strongly right-skewed. But in this heart,
`pctMito` runs from near 0 to 99%. On a log scale the cells near zero
stretch the spread so far that 3 MADs above the median lands past 100%.

On the original scale the cutoff is usable:

``` r

flagged <- flagQC(flagged, metrics = c(pctMito = "higher"), log = FALSE)
#> ℹ Flagging outliers in 1 metric at 3 MADs from the median, across all cells.
#> ℹ pctMito: 149/1011 cells flagged (14.74%).
#> • upper: > 42.29 -- 149 cells
#> ✔ Added 2 outlier columns; summary in `flagged@misc$flagQC`.
```

Check the thresholds in the console output before you trust the flags. A
cutoff outside the range a metric can take is the clearest sign the
scale is wrong for it.

## Estimating thresholds per sample

This is the argument that matters most. By default
[`flagQC()`](https://mikemartinez99.github.io/rscprep/reference/flagQC.md)
estimates one threshold across every cell in the object. If your samples
differ in sequencing depth – and they usually do – that single threshold
describes none of them.

The three pseudo-samples have median depths of roughly 700, 2,340, and
5,350 UMIs. Globally:

``` r

global <- flagQC(seurat, c(nFeature_RNA = "lower"))
#> ℹ Flagging outliers in 1 metric at 3 MADs from the median (log2 scale), across all cells.
#> ℹ nFeature_RNA: 14/1011 cells flagged (1.38%).
#> • lower: < 10.91 -- 14 cells
#> ✔ Added 1 outlier column; summary in `seurat@misc$flagQC`.

table(global$Ident, flagged = global$nFeature_RNA_3_lower_outliers)
#>           flagged
#>            FALSE TRUE
#>   Sample_A   343   13
#>   Sample_B   340    1
#>   Sample_C   314    0
```

One cutoff, set mostly by the middle of the pooled data. It sits low for
the deep sample and high for the shallow one, so its flags land on
`Sample_A` because that sample was sequenced shallowly, not because its
cells are worse. The deep sample’s poorest cells clear it easily.

Now per sample:

``` r

seurat <- flagQC(seurat, c(nFeature_RNA = "lower"), group_by = "Ident")
#> ℹ Flagging outliers in 1 metric at 3 MADs from the median (log2 scale), within 3 groups of Ident.
#> ℹ nFeature_RNA: 46/1011 cells flagged (4.55%) across 3 groups.
#> • lower: < 3.781-87.11 -- 46 cells
#> ✔ Added 1 outlier column; summary in `seurat@misc$flagQC`.

table(seurat$Ident, flagged = seurat$nFeature_RNA_3_lower_outliers)
#>           flagged
#>            FALSE TRUE
#>   Sample_A   353    3
#>   Sample_B   333    8
#>   Sample_C   279   35
```

Three thresholds instead of one, each set against its own sample’s
distribution. The flags move to the cells that are poor *for their
sample*, which in this case means mostly the deep one. The console
prints the *range* across groups rather than a line per group, so a
twenty-sample object stays readable; the per-group numbers are in the
summary table below.

`"Ident"` is the sample name
[`seuratFrom10x()`](https://mikemartinez99.github.io/rscprep/reference/seuratFrom10x.md)
sets, but any discrete metadata column works: a capture batch, a 10x
lane, a tissue. Group on the technical variable you want the threshold
to absorb.

`pctMito` is a ratio, so it does not scale with depth the way
`nFeature_RNA` does, and depth is not a reason to group it. Flag it
across all cells:

``` r

seurat <- flagQC(seurat, c(pctMito = "higher"), log = FALSE)
#> ℹ Flagging outliers in 1 metric at 3 MADs from the median, across all cells.
#> ℹ pctMito: 149/1011 cells flagged (14.74%).
#> • upper: > 42.29 -- 149 cells
#> ✔ Added 2 outlier columns; summary in `seurat@misc$flagQC`.
```

## Reading the summary

Every call records what it did.
[`getOutlierSummary()`](https://mikemartinez99.github.io/rscprep/reference/getOutlierSummary.md)
reads it back:

``` r

qc <- getOutlierSummary(seurat)
qc[, c("metric", "group", "lower", "upper", "n_flagged", "pct_flagged")]
#>         metric    group     lower    upper n_flagged pct_flagged
#> 1 nFeature_RNA Sample_A  3.780666       NA         3   0.8426966
#> 2 nFeature_RNA Sample_B 20.495647       NA         8   2.3460411
#> 3 nFeature_RNA Sample_C 87.106975       NA        35  11.1464968
#> 4      pctMito     <NA>        NA 42.29134       149  14.7378833
```

One row per metric per group. Two things to notice:

- Thresholds are on the **metric’s own scale**, even though the
  `nFeature_RNA` estimate was made on a log scale. They can be compared
  against the column directly.
- The unenforced tail is `NA`, not an unapplied number. `nFeature_RNA`
  was flagged `"lower"` only, so its `upper` is blank rather than a
  cutoff nothing was tested against. The same applies to `n_lower` /
  `n_upper`.

Which group is losing the most cells?

``` r

qc[order(-qc$pct_flagged), c("metric", "group", "n_flagged", "pct_flagged")]
#>         metric    group n_flagged pct_flagged
#> 4      pctMito     <NA>       149  14.7378833
#> 3 nFeature_RNA Sample_C        35  11.1464968
#> 2 nFeature_RNA Sample_B         8   2.3460411
#> 1 nFeature_RNA Sample_A         3   0.8426966
```

The two calls above accumulated into this one table: `nFeature_RNA` per
sample, then `pctMito` across all cells. It is a record written at flag
time, not a live view – cell counts describe the object as it stood
then, which is the point, since being able to say afterwards how many
cells were dropped and why is the reason to keep it.

## Groups too small to estimate from

Per-group estimation creates a failure mode global estimation does not
have. A group of two cells has no usable spread, and a group where over
half the cells share one value has a MAD of zero. Both are skipped
rather than given a meaningless threshold. Here, a made-up `Lane` column
puts two cells in a lane of their own:

``` r

lanes <- seurat
lanes$Lane <- ifelse(seq_len(ncol(lanes)) <= 2, "Lane_3",
                     rep(c("Lane_1", "Lane_2"), length.out = ncol(lanes)))

lanes <- flagQC(lanes, c(nCount_RNA = "lower"), group_by = "Lane")
#> ℹ Flagging outliers in 1 metric at 3 MADs from the median (log2 scale), within 3 groups of Lane.
#> Warning: ! nCount_RNA: no threshold for 1 of 3 Lane groups.
#> ℹ Skipped: Lane_3 (fewer than 3 cells)
#> ✖ Cells in that group are left as NA, not FALSE.
#> ℹ nCount_RNA: 0/1009 cells flagged (0%) across 2 groups.
#> • lower: < 16.33-18.42 -- 0 cells
#> ✔ Added 3 outlier columns; summary in `lanes@misc$flagQC`.
```

One warning naming every skipped group, not one per group. Skipped
groups still get a row in the summary table, with `NA` thresholds and an
`NA` `column`, so nothing silently disappears.

That `NA`-not-`FALSE` distinction runs throughout. A cell is left `NA`
when it was never judged – missing metric, missing `group_by` value, or
a skipped group – so an unjudged cell is never mistaken for one that
passed. It matters at the filtering step below.

## Actually filtering

[`flagQC()`](https://mikemartinez99.github.io/rscprep/reference/flagQC.md)
writes logical columns and stops. Filtering is yours:

``` r

flags <- grep("_outliers$", colnames(seurat@meta.data), value = TRUE)
flags
#> [1] "nFeature_RNA_3_lower_outliers" "pctMito_3_higher_outliers"

# `%in% TRUE` treats an unjudged NA as "not an outlier"; use `isTRUE`-style
# logic deliberately rather than letting NA propagate through `|`
keep     <- !Reduce(`|`, lapply(seurat@meta.data[flags], function(x) x %in% TRUE))
filtered <- seurat[, which(keep)]

table(seurat$Ident, keep)
#>           keep
#>            FALSE TRUE
#>   Sample_A    58  298
#>   Sample_B    59  282
#>   Sample_C    74  240
```

Cross-tabulate before you subset. A sample losing far more than its
neighbours is worth understanding rather than filtering: it may have
failed, or it may be a genuinely different tissue where the thresholds
do not transfer.

You do not have to trust every flag. The columns are independent, so
filter on the subset you believe:

``` r

filtered <- seurat[, which(!seurat$pctMito_3_higher_outliers %in% TRUE)]
ncol(filtered)
#> [1] 862
```

## Tuning

``` r

# more aggressive
seurat <- flagQC(seurat, c(pctMito = "higher"), nmads = 2.5, log = FALSE)

# keep all but the extremes
seurat <- flagQC(seurat, c(pctMito = "higher"), nmads = 5, log = FALSE)

# a metric that is already symmetric does not want the log transform
seurat <- flagQC(seurat, c(log10GenesPerUMI = "lower"), log = FALSE)
```

``` r

colSums(seurat@meta.data[grep("^pctMito_", colnames(seurat@meta.data))],
        na.rm = TRUE)
#>   pctMito_3_higher_outliers pctMito_2.5_higher_outliers 
#>                         149                         186 
#>   pctMito_5_higher_outliers 
#>                          45
```

`nmads` is in the column name, so `nmads = 2.5`, `3` and `5` produce
separate columns and can be compared side by side. `group_by` is **not**
in the name, so re-running one metric with different grouping overwrites
the previous flags – the summary table is updated to match, so it never
describes a column that is no longer there.

`log = TRUE` is the default and is right for count metrics, which are
right-skewed enough that the median and MAD are misleading on the raw
scale. Turn it off for anything already roughly symmetric, and for
bounded metrics like `pctMito` whose low end crowds against zero.

## Where this leaves you

``` r

seurat <- mergeSeurat(seuratFrom10x(manifest = manifest))
seurat <- diagnoseQC(seurat, organism = "mmu")
seurat <- flagQC(
  seurat,
  metrics  = c(nCount_RNA = "lower", nFeature_RNA = "lower"),
  group_by = "Ident"
)
seurat <- flagQC(seurat, c(pctMito = "higher"), log = FALSE)

getOutlierSummary(seurat)          # what would be cut, and from where
write.csv(getOutlierSummary(seurat), "qc/outlier_summary.csv", row.names = FALSE)
```

Nothing has been removed at this point. The metrics are on the object,
the flags are on the object, and the record of how the thresholds were
chosen is on the object. Deciding what to drop is the next step, and it
is a judgement call – which is why these functions stop short of making
it.
