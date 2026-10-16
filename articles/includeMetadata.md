# Assigning samples to treatment groups

``` r

library(rscprep)
library(Seurat)
```

Treatment group, genotype, timepoint – the things you actually want to
compare – have to get attached to the objects somewhere.
[`includeMetadata()`](https://mikemartinez99.github.io/rscprep/reference/includeMetadata.md)
does that from a sample-level table, keyed on the names of the list
[`seuratFrom10x()`](https://mikemartinez99.github.io/rscprep/reference/seuratFrom10x.md)
returns.

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

The treatment labels below are invented for the example. All three
pseudo-samples are the same heart.

## Why not just put it in the manifest?

Because the manifest is a different kind of table. Its rows are Cell
Ranger **runs** – one per `outs/` directory. Treatment is a property of
a **sample**. For a plain `cellranger count` run those happen to
coincide, so it is tempting to add a `Treatment` column and be done:

``` r

manifest <- data.frame(
  Path        = "/path/to/Sample_1/outs",
  multi = FALSE,
  Treatment   = "Control"     # works here, and only here
)
```

But a multiplexed run breaks it. One row of the manifest, one `outs/`
directory, and eight demuxed samples inside it – possibly four controls
and four treated. There is no way to write four values in one cell that
isn’t a delimited string you then have to parse and validate.

So the two tables stay separate, at their own grain:

- **manifest** – run-level: `Path`, `multi`, and legitimately things
  like `Batch` or `RunDate` that really do apply to the whole run.
- **sample sheet** – sample-level: one row per sample, however many runs
  they came from.

By the time
[`seuratFrom10x()`](https://mikemartinez99.github.io/rscprep/reference/seuratFrom10x.md)
hands you a list, every entry is one sample and the
multiplexed/non-multiplexed distinction has been erased. That is the
right moment to attach sample-level metadata, and it is the same code
path either way.

## The workflow

Start from the list and read its names:

``` r

names(seurat_list)
#> [1] "Sample_A" "Sample_B" "Sample_C"
```

For a `cellranger multi` run those are the folder names under
`per_sample_outs/`. Those names are what you key on. Write a sheet
against them:

``` r

sheet <- data.frame(
  Sample    = c("Sample_A", "Sample_B", "Sample_C"),
  Treatment = c("Control", "Bleomycin", "Bleomycin"),
  Timepoint = c("D21", "D21", "D21")
)

annotated <- includeMetadata(seurat_list, sheet)
#> ✔ Added 2 metadata columns to 3 objects: "Treatment" and "Timepoint"
```

In practice the sheet lives in a CSV next to your analysis, filled in by
whoever ran the experiment and kept under version control:

``` r

sheet     <- read.csv("samplesheet.csv")
annotated <- includeMetadata(seurat_list, sheet)
```

Each column is broadcast across every cell of its sample, so after
merging you can group by it directly:

``` r

merged <- mergeSeurat(annotated)
#> ℹ Merging 3 seurat objects...
#> ✔ Successfully merged 3 seurat objects.
#> An object of class Seurat 
#> 2092 features across 1011 samples within 1 assay 
#> Active assay: RNA (2092 features, 0 variable features)
#>  1 layer present: counts

table(merged$Ident, merged$Treatment)
#>           
#>            Bleomycin Control
#>   Sample_A         0     356
#>   Sample_B       341       0
#>   Sample_C       314       0
```

You don’t have to know the sample names in advance – `per_sample_outs/`
does. Load first, read `names(seurat_list)`, and write the sheet against
what is actually there. Discovery is automatic; the assignment stays
explicit.

### Identifying samples by rownames

A `Sample` column is the primary key because it survives a CSV round
trip. Rownames work too, and are used whenever there is no `Sample`
column:

``` r

sheet <- data.frame(Treatment = c("Control", "Bleomycin", "Bleomycin"))
rownames(sheet) <- c("Sample_A", "Sample_B", "Sample_C")

annotated <- includeMetadata(seurat_list, sheet)
#> ✔ Added 1 metadata column to 3 objects: "Treatment"
```

If both are present the column wins. A data frame with neither – so,
default `1, 2, 3` rownames – is rejected rather than matched
positionally, because lining rows up by position against a list whose
order you did not choose is a silent way to mislabel an entire
experiment.

## Every column, every sample

The rule
[`includeMetadata()`](https://mikemartinez99.github.io/rscprep/reference/includeMetadata.md)
enforces is that a column, once present, is populated for **all**
samples. A sample with no row is an error. So is an `NA`, an empty cell,
or a whitespace-only cell:

``` r

sheet <- data.frame(
  Sample    = c("Sample_A", "Sample_B", "Sample_C"),
  Treatment = c("Control", NA, "Bleomycin"),
  Timepoint = c(NA, "D21", "")
)

includeMetadata(seurat_list, sheet)
#> Error in `includeMetadata()`:
#> ! Every column of `metadata` must be populated for every sample.
#> ✖ Treatment: missing for "Sample_B"
#> ✖ Timepoint: missing for "Sample_A" and "Sample_C"
#> ℹ A column left blank for some samples becomes `NA` for their cells after
#>   `mergeSeurat()`.
```

Every offending column is reported at once, so a forty-sample sheet gets
fixed in one pass rather than forty.

This strictness exists because the alternative is silent.
`Seurat::merge()` will happily combine an object carrying a `Treatment`
column with one that doesn’t – it fills the gap with `NA` and says
nothing:

``` r

annotated   <- seurat_list$Sample_A
unannotated <- seurat_list$Sample_B
annotated$Treatment <- "Control"

merged <- merge(annotated, unannotated, add.cell.id = c("Sample_A", "Sample_B"))

table(merged$Treatment, useNA = "ifany")
#> 
#> Control    <NA> 
#>     356     341
```

Those `NA` cells then quietly drop out of every grouped comparison
downstream. By the time you notice a group is short, the cause is
several steps back. Requiring completeness up front makes the state
unreachable.

## Master sheets

Rows that match no loaded object are ignored with a warning, so one
sheet can describe a whole experiment while you work on a subset:

``` r

# the sheet describes five samples; only three are loaded
master_sheet <- data.frame(
  Sample    = c("Sample_A", "Sample_B", "Sample_C", "Sample_D", "Sample_E"),
  Treatment = c("Control", "Bleomycin", "Bleomycin", "Control", "Bleomycin")
)

annotated <- includeMetadata(seurat_list, master_sheet)
#> Warning: Ignoring 2 rows of `metadata` matching no object in `objList`: "Sample_D" and
#> "Sample_E"
#> ✔ Added 1 metadata column to 3 objects: "Treatment"
```

A typo is still caught, because a mistyped name in the sheet leaves some
object without a row – and that is a hard error, not a warning.

## Two smaller rules

**Existing columns are protected.** Re-running against an edited sheet
errors unless you opt in:

``` r

corrected_sheet <- master_sheet
corrected_sheet$Treatment[2] <- "Control"

includeMetadata(annotated, corrected_sheet)
#> Warning: Ignoring 2 rows of `metadata` matching no object in `objList`: "Sample_D" and
#> "Sample_E"
#> Error in `includeMetadata()`:
#> ! 1 column in `metadata` is already present in `objList`: "Treatment"
#> ℹ Set `overwrite = TRUE` to replace them.

annotated <- includeMetadata(annotated, corrected_sheet, overwrite = TRUE)
#> Warning: Ignoring 2 rows of `metadata` matching no object in `objList`: "Sample_D" and
#> "Sample_E"
#> ✔ Added 1 metadata column to 3 objects: "Treatment"
```

`Ident` is reserved and cannot be set at all, even with
`overwrite = TRUE` – it identifies the sample of origin and is what
[`mergeSeurat()`](https://mikemartinez99.github.io/rscprep/reference/mergeSeurat.md)
verifies after merging.

**Column names must be syntactically valid.** Seurat quietly runs
[`make.names()`](https://rdrr.io/r/base/make.names.html) on metadata
column names, so `Cell Type` would land as `Cell.Type` and stop matching
the header in your sheet.
[`includeMetadata()`](https://mikemartinez99.github.io/rscprep/reference/includeMetadata.md)
refuses instead, and tells you what the name would have become. Note
that [`read.csv()`](https://rdrr.io/r/utils/read.table.html) already
does this conversion for you by default, so it rarely comes up.

## A note on factors

Column types are preserved, factors included – each object keeps the
sheet’s full level set, not just the level it happens to carry, so
per-sample subsets stay comparable.

That survives on the list but **not** through a merge: `Seurat::merge()`
flattens factors to character. If you want an ordered factor on the
merged object, re-apply it afterwards:

``` r

merged <- mergeSeurat(annotated)
#> ℹ Merging 3 seurat objects...
#> ✔ Successfully merged 3 seurat objects.
#> An object of class Seurat 
#> 2092 features across 1011 samples within 1 assay 
#> Active assay: RNA (2092 features, 0 variable features)
#>  1 layer present: counts
merged$Treatment <- factor(merged$Treatment, levels = c("Control", "Bleomycin"))
```
