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

``` r

seurat_list <- seuratFrom10x(
  path        = "/path/to/Pool1/outs",
  multi = TRUE
)

names(seurat_list)
#> [1] "Bleo_18_Slim_10" "Bleo_22_Slim_14" "Ctrl_04_Slim_02" "Ctrl_09_Slim_07"
```

Those names are what you key on. Write a sheet against them:

``` r

sheet <- data.frame(
  Sample    = c("Bleo_18_Slim_10", "Bleo_22_Slim_14",
                "Ctrl_04_Slim_02", "Ctrl_09_Slim_07"),
  Treatment = c("Bleomycin", "Bleomycin", "Control", "Control"),
  Timepoint = c("D21", "D21", "D21", "D21")
)

annotated <- includeMetadata(seurat_list, sheet)
#> ✔ Added 2 metadata columns to 4 objects: "Treatment" and "Timepoint"
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

table(merged$Ident, merged$Treatment)
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

sheet <- data.frame(Treatment = c("Control", "Bleomycin"))
rownames(sheet) <- c("Sample_A", "Sample_B")

annotated <- includeMetadata(seurat_list, sheet)
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
#> Error:
#> Every column of `metadata` must be populated for every sample.
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
#> Control    <NA>
#>    4821    5104
```

Those `NA` cells then quietly drop out of every grouped comparison
downstream. By the time you notice a group is short, the cause is
several steps back. Requiring completeness up front makes the state
unreachable.

## Master sheets

Rows that match no loaded object are ignored with a warning, so one
sheet can describe a whole experiment while you work on a subset:

``` r

# sheet covers all 40 samples; only Pool1's 4 are loaded
annotated <- includeMetadata(seurat_list, master_sheet)
#> ! Ignoring 36 rows of `metadata` matching no object in `objList`: ...
#> ✔ Added 2 metadata columns to 4 objects: "Treatment" and "Timepoint"
```

A typo is still caught, because a mistyped name in the sheet leaves some
object without a row – and that is a hard error, not a warning.

## Two smaller rules

**Existing columns are protected.** Re-running against an edited sheet
errors unless you opt in:

``` r

includeMetadata(annotated, corrected_sheet, overwrite = TRUE)
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
merged$Treatment <- factor(merged$Treatment, levels = c("Control", "Bleomycin"))
```
