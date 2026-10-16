# Merging per-sample Seurat objects

``` r

library(rscprep)
library(Seurat)
```

[`seuratFrom10x()`](https://mikemartinez99.github.io/rscprep/reference/seuratFrom10x.md)
hands back a named list – one Seurat object per sample.
[`mergeSeurat()`](https://mikemartinez99.github.io/rscprep/reference/mergeSeurat.md)
is the other half of that pairing: it collapses the list into a single
object you can take through a normalisation and integration workflow,
without losing track of which cell came from which sample.

## The basic pipeline

``` r

seurat_list <- seuratFrom10x(
  path        = "/path/to/Pool1/outs",
  multi = TRUE
)

merged <- mergeSeurat(seurat_list)
```

That is the whole intended usage: the list
[`seuratFrom10x()`](https://mikemartinez99.github.io/rscprep/reference/seuratFrom10x.md)
returns is already in the shape
[`mergeSeurat()`](https://mikemartinez99.github.io/rscprep/reference/mergeSeurat.md)
expects, so nothing needs renaming or reordering in between.

## What the list names are for

The names of `objList` are not decoration – they are passed to
[`merge()`](https://rdrr.io/r/base/merge.html) as `add.cell.id`, so each
becomes a prefix on that sample’s cell barcodes:

``` r

colnames(seurat_list$Sample_A)[1:2]
#> [1] "AAACCCAAGACCTTTG-1" "AAACCCACAGCTATTG-1"

colnames(merged)[1:2]
#> [1] "Sample_A_AAACCCAAGACCTTTG-1" "Sample_A_AAACCCACAGCTATTG-1"
```

This matters because 10x barcodes are drawn from a fixed whitelist, so
the same barcode routinely appears in several samples. Prefixing keeps
them distinct and keeps every cell traceable back to its origin. It is
also why names must be present and unique –
[`mergeSeurat()`](https://mikemartinez99.github.io/rscprep/reference/mergeSeurat.md)
refuses to merge a list that is unnamed, has blank or `NA` names, or
repeats a name.

The `Ident` column set by
[`seuratFrom10x()`](https://mikemartinez99.github.io/rscprep/reference/seuratFrom10x.md)
is carried through the merge, so you can still split or group by sample
afterwards:

``` r

table(merged$Ident)
#> Sample_A Sample_B
#>     4821     5104
```

## What gets checked

The list is validated before anything is merged, so a bad batch fails
with a message naming the offending sample rather than producing a
subtly wrong object.
[`mergeSeurat()`](https://mikemartinez99.github.io/rscprep/reference/mergeSeurat.md)
requires that `objList`:

- is non-empty;
- contains only structurally valid Seurat objects;
- is fully named, with no blank, `NA` or duplicated names;
- contains no zero-cell objects;
- agrees on
  [`DefaultAssay()`](https://satijalab.github.io/seurat-object/reference/DefaultAssay.html)
  across every element.

The assay check is about *agreement*, not about `RNA` specifically – a
list of `ADT` objects merges fine, a list mixing `RNA` and `ADT` does
not.

Two more checks run after the merge, to catch a merge that silently lost
data: every sample must still be present in the merged `Ident` column,
and the merged cell count must equal the sum of the inputs’.

## Two behaviours worth knowing

A single-element list short-circuits. The lone object is returned as-is,
with a message – unmerged, and with its barcodes left unprefixed:

``` r

one <- seuratFrom10x(path = "/path/to/Sample_1/outs", multi = FALSE)

merged <- mergeSeurat(one)
#> i 'one' has only one object; returning it unmerged.
```

And the merged object is returned **invisibly**, after being printed. So
a bare call shows you the result without flooding the console twice, but
you still need to assign it to keep it:

``` r

mergeSeurat(seurat_list)      # prints a summary; result is not kept
merged <- mergeSeurat(seurat_list)   # assign to keep it
```
