# Attach Sample-Level Metadata to a List of Seurat Objects

Joins a sample-level table – treatment group, genotype, timepoint, and
so on – onto a named list of
[Seurat](https://satijalab.org/seurat/reference/reexports.html) objects,
typically the output of
[`seuratFrom10x()`](https://mikemartinez99.github.io/rscprep/reference/seuratFrom10x.md).
Each row of `metadata` describes one sample; every other column becomes
a metadata column on that sample's object, broadcast across all of its
cells.

## Usage

``` r
includeMetadata(objList, metadata, overwrite = FALSE)
```

## Arguments

- objList:

  Named list of
  [Seurat](https://satijalab.org/seurat/reference/reexports.html)
  objects, one per sample, as returned by
  [`seuratFrom10x()`](https://mikemartinez99.github.io/rscprep/reference/seuratFrom10x.md).

- metadata:

  Data frame with one row per sample. Sample names come from a `Sample`
  column if present, otherwise from the rownames, and must match
  `names(objList)`. Every remaining column is attached.

- overwrite:

  Logical. Whether to replace metadata columns that already exist on the
  objects. Defaults to `FALSE`, which errors instead. `Ident` is
  reserved and cannot be set either way – it identifies the sample of
  origin and is what
  [`mergeSeurat()`](https://mikemartinez99.github.io/rscprep/reference/mergeSeurat.md)
  verifies after merging.

## Value

`objList` with the new columns attached, in the same order and with the
same names. Column types are preserved, including factor levels – though
note that `Seurat::merge()` flattens factors to character, so a factor
column survives on the list but not through
[`mergeSeurat()`](https://mikemartinez99.github.io/rscprep/reference/mergeSeurat.md).

## Details

This is deliberately separate from the manifest
[`seuratFrom10x()`](https://mikemartinez99.github.io/rscprep/reference/seuratFrom10x.md)
accepts. That manifest is a *run*-level table – one row per Cell Ranger
`outs/` directory – and a multiplexed run expands into many samples, so
a single manifest row cannot carry per-sample assignments. `metadata` is
a *sample*-level table instead, keyed on the names of `objList`. That
works identically for multiplexed and non-multiplexed input, because by
the time the list exists the distinction is already gone.

Samples are identified by a `Sample` column if `metadata` has one, and
by its rownames otherwise. A `Sample` column takes precedence, and is
used as the key rather than attached as metadata.

Every column must be populated for *every* object in `objList`: a sample
with no row in `metadata`, or an `NA` or blank value in any column, is
an error. This is what keeps the result safe to
[`mergeSeurat()`](https://mikemartinez99.github.io/rscprep/reference/mergeSeurat.md).
[`merge()`](https://rdrr.io/r/base/merge.html) will happily combine an
object carrying a `Treatment` column with one that lacks it, silently
filling the gap with `NA` for every cell of the un-annotated sample;
requiring completeness up front makes that state unreachable.

Rows matching no object are ignored with a warning, so one master sheet
can cover a whole experiment while you work on a subset of it. A
mistyped sample name is still caught, since it also leaves some object
without a row.

## See also

[`seuratFrom10x()`](https://mikemartinez99.github.io/rscprep/reference/seuratFrom10x.md),
which produces the named list this expects, and
[`mergeSeurat()`](https://mikemartinez99.github.io/rscprep/reference/mergeSeurat.md),
which consumes the annotated result.

## Examples

``` r
if (FALSE) { # \dontrun{
seurat_list <- seuratFrom10x(
  path        = "/path/to/Pool1/outs",
  multi = TRUE
)

# one row per demuxed sample, however many the run produced
sheet <- data.frame(
  Sample    = c("Sample_A", "Sample_B"),
  Treatment = c("Control", "Bleomycin"),
  Timepoint = c("D0", "D21")
)

annotated <- includeMetadata(seurat_list, sheet)
merged    <- mergeSeurat(annotated)

table(merged$Ident, merged$Treatment)
} # }
```
