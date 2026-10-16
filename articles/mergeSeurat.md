# Merging per-sample Seurat objects

``` r

library(rscprep)
library(Seurat)
```

[`seuratFrom10x()`](https://mikemartinez99.github.io/rscprep/reference/seuratFrom10x.md)
hands back a named list: one Seurat object per sample.
[`mergeSeurat()`](https://mikemartinez99.github.io/rscprep/reference/mergeSeurat.md)
is the other half of that pairing: it collapses the list into a single
merged object you can take through a normalisation and integration
workflow, without losing track of which cell came from which sample.

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

The names of `objList` are not decoration. They are passed to
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
also why names must be present and unique.
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

## Seurat 4 and Seurat 5

[`mergeSeurat()`](https://mikemartinez99.github.io/rscprep/reference/mergeSeurat.md)
works under both versions and detects which one is installed. The
difference it handles is what
[`merge()`](https://rdrr.io/r/base/merge.html) returns. Under Seurat 4,
every sample’s counts are bound into one matrix. Under Seurat 5, each
sample keeps its own layer: `counts.Sample_A`, `counts.Sample_B`, and so
on.

By default
[`mergeSeurat()`](https://mikemartinez99.github.io/rscprep/reference/mergeSeurat.md)
joins those layers, so the merged object has the same shape whichever
version made it:

``` r

merged <- mergeSeurat(seurat_list)
merged
#> An object of class Seurat
#> 33538 features across 9925 samples within 1 assay
#> Active assay: RNA (33538 features, 0 variable features)
#>  1 layer present: counts
```

Seurat 5’s own integration workflow,
[`IntegrateLayers()`](https://satijalab.org/seurat/reference/IntegrateLayers.html),
needs the layers split by sample. Ask for that with
`joinLayers = FALSE`:

``` r

merged_split <- mergeSeurat(seurat_list, joinLayers = FALSE)
merged_split
#> An object of class Seurat
#> 33538 features across 9925 samples within 1 assay
#> Active assay: RNA (33538 features, 0 variable features)
#>  2 layers present: counts.Sample_A, counts.Sample_B
```

What happens in each case:

| Situation | Result |
|----|----|
| Seurat 5, `Assay5` assays, `joinLayers = TRUE` (default) | Layers joined into one `counts` layer |
| Seurat 5, `joinLayers = FALSE` | Layers left split, one per sample |
| Seurat 4 | Nothing to join; [`merge()`](https://rdrr.io/r/base/merge.html) never splits, and `joinLayers` is ignored |
| A list of one object | Returned unmerged, so there are no layers to join |

Joining is decided assay by assay, from each assay’s class rather than
from the installed version alone. An object with several assays has only
its `Assay5` assays joined.
