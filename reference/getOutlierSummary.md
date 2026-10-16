# Retrieve the Outlier Summary Table from an Object

Returns the table
[`flagQC()`](https://mikemartinez99.github.io/rscprep/reference/flagQC.md)
recorded when it flagged the object: one row per metric per group,
holding the thresholds that were applied and how many cells each one
caught. It is a plain data frame, so it can be sorted, subset, written
to CSV or joined against a sample sheet like any other.

## Usage

``` r
getOutlierSummary(obj)
```

## Arguments

- obj:

  A [Seurat](https://satijalab.org/seurat/reference/reexports.html)
  object that has been through
  [`flagQC()`](https://mikemartinez99.github.io/rscprep/reference/flagQC.md).

## Value

A data frame with one row per metric per group. See the *Value* section
of
[`flagQC()`](https://mikemartinez99.github.io/rscprep/reference/flagQC.md)
for what each column means. Aborts if the object has never been flagged,
rather than returning `NULL` for a caller to trip over.

## Details

The table is a record of what was flagged, not a live view of the
object. It is written when
[`flagQC()`](https://mikemartinez99.github.io/rscprep/reference/flagQC.md)
runs and is not recalculated afterwards, so cell counts describe the
object as it stood at that moment. Subsetting the object does not update
them – which is the point, since the whole reason to keep the table is
to be able to say afterwards how many cells were dropped and why.

Repeated
[`flagQC()`](https://mikemartinez99.github.io/rscprep/reference/flagQC.md)
calls accumulate into one table, so an object flagged on counts globally
and then on mitochondrial content per sample yields a single table
describing both. Rows for a flag column that was overwritten are
replaced rather than duplicated, so every row corresponds to a metadata
column that is actually present.

## See also

[`flagQC()`](https://mikemartinez99.github.io/rscprep/reference/flagQC.md),
which writes this table, and
[`diagnoseQC()`](https://mikemartinez99.github.io/rscprep/reference/diagnoseQC.md),
which calculates the metrics it describes.

## Examples

``` r
if (FALSE) { # \dontrun{
seurat <- diagnoseQC(seurat, organism = "mmu")
seurat <- flagQC(
    seurat,
    metrics  = c(nCount_RNA = "lower", pctMito = "higher"),
    group_by = "Ident"
)

qc <- getOutlierSummary(seurat)
qc

# the thresholds that were actually applied, per sample
qc[, c("metric", "group", "lower", "upper", "n_lower", "n_upper")]

# which samples lose the most cells, and to which metric?
qc[order(-qc$pct_flagged), c("metric", "group", "n_flagged", "pct_flagged")]

# groups no threshold could be estimated for
qc[is.na(qc$column), c("metric", "group", "n_cells")]

# keep it alongside the filtered object
write.csv(qc, "qc/outlier_summary.csv", row.names = FALSE)
} # }
```
