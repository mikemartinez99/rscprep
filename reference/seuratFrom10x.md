# Read 10x Genomics Gene Expression Sample(s) into Seurat Objects

Reads one or more samples' filtered feature barcode matrix H5 files from
10x Genomics Cell Ranger output and creates one
[Seurat](https://satijalab.org/seurat/reference/reexports.html) object
per sample. Supports both `cellranger count` output (one sample per run)
and `cellranger multi` output (one or more demuxed samples under
`outs/per_sample_outs/`), specified either as a single path or as a
manifest describing multiple runs.

## Usage

``` r
seuratFrom10x(
  manifest = NULL,
  path = NULL,
  multi = NULL,
  min.cells = 10,
  min.features = 50,
  multiplexed = NULL
)
```

## Arguments

- manifest:

  Data frame with columns `Path` (character; path to each run's `outs/`
  directory) and `multi` (logical). Rows may mix the two pipelines
  freely. Mutually exclusive with `path`.

- path:

  Character string. Path to a single Cell Ranger `outs/` directory.
  Mutually exclusive with `manifest`.

- multi:

  Logical. Whether `path` was produced by `cellranger multi` rather than
  `cellranger count` – see the section above. Required when `path` is
  supplied; ignored when `manifest` is supplied, since the pipeline is
  read per-row from `manifest$multi` instead.

- min.cells:

  Integer. Include features detected in at least this many cells. Passed
  to
  [`Seurat::CreateSeuratObject()`](https://satijalab.org/seurat/reference/reexports.html).

- min.features:

  Integer. Include cells where at least this many features are detected.
  Passed to
  [`Seurat::CreateSeuratObject()`](https://satijalab.org/seurat/reference/reexports.html).

- multiplexed:

  Deprecated. The former name of `multi`; still accepted, with a
  warning.

## Value

A named list of
[Seurat](https://satijalab.org/seurat/reference/reexports.html) objects,
one per sample, each with the `RNA` assay and an `Ident` metadata column
set to the sample name. List names match sample identifiers. Samples
whose expected H5 file is missing are dropped from the manifest-mode
result with a warning; in `path` mode, a missing file raises an error
instead.

## Details

Exactly one of `manifest` or `path` must be supplied.

## Which pipeline wrote the run

`multi` selects the output *layout*, and the two Cell Ranger pipelines
produce different ones:

- `multi = FALSE`:

  `cellranger count`, which writes `outs/filtered_feature_bc_matrix.h5`
  for its single sample.

- `multi = TRUE`:

  `cellranger multi`, which writes
  `outs/per_sample_outs/<sample>/count/sample_filtered_feature_bc_matrix.h5`
  for each sample it demuxed.

The distinction is which pipeline ran, not how many samples came out of
it. `cellranger multi` writes `per_sample_outs/` even when it
demultiplexed a single sample – the usual case for Flex and 5' runs – so
those need `multi = TRUE` despite not being multiplexed in any
biological sense. That is why the argument is no longer called
`multiplexed`.

## Examples

``` r
if (FALSE) { # \dontrun{
# A `cellranger count` run
seurat_list <- seuratFrom10x(
  path  = "/path/to/Sample_1/outs",
  multi = FALSE
)

# A `cellranger multi` run, expanding to one object per demuxed sample
seurat_list <- seuratFrom10x(
  path  = "/path/to/Pool1/outs",
  multi = TRUE
)

# Manifest describing several runs from either pipeline
manifest <- data.frame(
  Path  = c("/path/to/Sample_1/outs", "/path/to/Pool1/outs"),
  multi = c(FALSE, TRUE)
)
seurat_list <- seuratFrom10x(manifest = manifest, min.cells = 3, min.features = 200)
} # }
```
