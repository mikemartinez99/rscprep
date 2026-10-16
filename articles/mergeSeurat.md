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

## The basic pipeline

With real data, the list comes straight from
[`seuratFrom10x()`](https://mikemartinez99.github.io/rscprep/reference/seuratFrom10x.md):

``` r

seurat_list <- seuratFrom10x(
  path  = "/path/to/Pool1/outs",
  multi = TRUE
)
```

Here it is the pseudo-sample list built above, and the merge is the same
either way:

``` r

merged <- mergeSeurat(seurat_list)
#> ℹ Merging 3 seurat objects...
#> ✔ Successfully merged 3 seurat objects.
#> An object of class Seurat 
#> 2092 features across 1011 samples within 1 assay 
#> Active assay: RNA (2092 features, 0 variable features)
#>  1 layer present: counts
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
#> [1] "AAACCCACACCAGTAT-1" "AAACCCAGTGGAACAC-1"

colnames(merged)[1:2]
#> [1] "Sample_A_AAACCCACACCAGTAT-1" "Sample_A_AAACCCAGTGGAACAC-1"
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
#> 
#> Sample_A Sample_B Sample_C 
#>      356      341      314
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

``` r

mergeSeurat(setNames(seurat_list, c("Sample_A", "Sample_A", "Sample_C")))
#> Error in `mergeSeurat()`:
#> ! 'setNames(seurat_list, c("Sample_A", "Sample_A", "Sample_C"))'
#>   contains duplicate names: Sample_A
```

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

one <- seuratFrom10x(path = outs, multi = FALSE)
#> ✔ "heart_1k_v3": 957 cells, 2089 features.
#> ✔ Loaded 1 sample.

merged_one <- mergeSeurat(one)
#> ℹ 'one' has only one object; returning it unmerged.
colnames(merged_one)[1:2]
#> [1] "AAACCCACACCAGTAT-1" "AAACCCAGTCACCTTC-1"
```

And the merged object is returned **invisibly**, after being printed. So
a bare call shows you the result without flooding the console twice, but
you still need to assign it to keep it:

``` r

mergeSeurat(seurat_list)             # prints a summary; result is not kept
merged <- mergeSeurat(seurat_list)   # assign to keep it
```

## Seurat 4 and Seurat 5

[`mergeSeurat()`](https://mikemartinez99.github.io/rscprep/reference/mergeSeurat.md)
detects which Seurat is installed. Under Seurat 4,
[`merge()`](https://rdrr.io/r/base/merge.html) binds every sample’s
counts into one matrix. Under Seurat 5 it keeps one layer per sample
instead (`counts.Sample_A`, `counts.Sample_B`, …).
[`mergeSeurat()`](https://mikemartinez99.github.io/rscprep/reference/mergeSeurat.md)
joins those layers by default, so the merged object has the same shape
under either version and code downstream does not need to know which one
made it.

Seurat 5’s own integration workflow,
[`IntegrateLayers()`](https://satijalab.org/seurat/reference/IntegrateLayers.html),
wants the layers split. Ask for that with `joinLayers = FALSE`; under
Seurat 4 the argument does nothing:

``` r

merged <- mergeSeurat(seurat_list, joinLayers = FALSE)
```
