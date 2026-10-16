# Reading 10x Genomics data into Seurat objects

``` r

library(rscprep)
library(Seurat)
```

[`seuratFrom10x()`](https://mikemartinez99.github.io/rscprep/reference/seuratFrom10x.md)
is a helper function to efficiently read in single cell datasets from
10X. This function supports both 10X CellRanger output formats: standard
single-sample `cellranger count` runs, and demuxed `cellranger multi`
runs.

You can point it at a single sample or directory, or hand it a manifest
file describing a batch of multiple runs from either pipeline or both.
It returns a named list of Seurat objects, one per sample, ready to
merge or process independently. Importantly, the package auto detects
which version of Seurat you have to account for `SeuratObject`
differences between v4 and v5.

The `multi` argument controls how the function reads in your data. It
tells the code which pipeline wrote the directory, `count` or `multi`.

## A `cellranger count` run

``` r

seurat_list <- seuratFrom10x(
  path  = "/path/to/Sample_1/outs",
  multi = FALSE
)

seurat_list$Sample_1
```

`path` should point to the Cell Ranger `outs/` directory. The function
then looks for `filtered_feature_bc_matrix.h5` directly inside it. The
sample name is inferred from the parent folder (`Sample_1`, one level up
from `outs`).

## A `cellranger multi` run

``` r

seurat_list <- seuratFrom10x(
  path  = "/path/to/Pool1/outs",
  multi = TRUE
)

names(seurat_list)
#> [1] "Bleo_18_Slim_10" "Bleo_22_Slim_14" ...
```

For `cellranger multi` output, one `outs/` directory holds every demuxed
sample under `per_sample_outs/<sample>/`.
[`seuratFrom10x()`](https://mikemartinez99.github.io/rscprep/reference/seuratFrom10x.md)
reads each sample’s `count/sample_filtered_feature_bc_matrix.h5`
separately and returns one Seurat object per demuxed sample. A single
function call here returns a list, not a single object.

## A batch of runs via manifest

``` r

manifest <- data.frame(
  Path = c(
    "/path/to/Sample_1/outs",
    "/path/to/Sample_2/outs",
    "/path/to/Pool1/outs"
  ),
  multi = c(FALSE, FALSE, TRUE)
)

seurat_list <- seuratFrom10x(manifest = manifest, min.cells = 3, min.features = 200)
```

The manifest needs a `Path` column (each row’s `outs/` directory) and a
logical `multi` column passed as a vector. Rows can mix the two
pipelines freely. If a sample’s expected H5 file is missing, that row is
skipped with a warning rather than failing the whole batch; single-path
mode fails immediately instead, since there’s nothing left to fall back
on.

## Choosing `multi`

`multi` picks between two directory layouts.

| `multi` | Pipeline | Where the filtered H5 lives |
|----|----|----|
| `FALSE` | `cellranger count` | `outs/filtered_feature_bc_matrix.h5` |
| `TRUE` | `cellranger multi` | `outs/per_sample_outs/<sample>/count/sample_filtered_feature_bc_matrix.h5` |

The question is which pipeline ran, not how many samples came out of it.
`cellranger multi` writes `per_sample_outs/` even when it demultiplexed
exactly one sample, so those runs need `multi = TRUE` despite not being
multiplexed in any biological sense. If you are unsure, look for
`per_sample_outs/` inside `outs/`: if it is there, the run needs
`multi = TRUE`.

**This argument used to be called `multiplexed`, which read as a claim
about the experiment rather than about the pipeline. The old name still
works and warns.**

## Notes

`min.cells` and `min.features` are passed straight through to
[`Seurat::CreateSeuratObject()`](https://satijalab.github.io/seurat-object/reference/CreateSeuratObject.html)
and apply uniformly across every sample in a batch.

The code auto-detects the Seurat version in your session and reads in as
either classic assay or Assay5 accordingly.
