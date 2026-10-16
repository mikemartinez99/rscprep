# Read 10x Genomics Raw (Unfiltered) Feature Barcode Matrices

The raw-matrix counterpart of
[`seuratFrom10x()`](https://mikemartinez99.github.io/rscprep/reference/seuratFrom10x.md).
Reads the *unfiltered* feature barcode matrix H5 – every barcode Cell
Ranger saw, not just the ones it called as cells – for one or more runs,
taking the same `manifest` / `path` / `multi` arguments.

## Usage

``` r
rawMatrixFrom10x(
  manifest = NULL,
  path = NULL,
  multi = NULL,
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

- multiplexed:

  Deprecated. The former name of `multi`; still accepted, with a
  warning.

## Value

A named list of sparse count matrices (features x barcodes), one per
run, named after the directory above each `outs/`. Runs whose raw H5 is
missing are dropped from the manifest-mode result with a warning; in
`path` mode a missing file raises an error instead.

## Details

This is the input a barcode rank plot needs. Cell calling is exactly the
step that separates the raw matrix from the filtered one, so a plot
meant to show *where that cutoff landed* cannot be drawn from the
filtered matrix: the empty droplets that form the left-hand cliff have
already been removed. Use
[`plotBarcodeRanks()`](https://mikemartinez99.github.io/rscprep/reference/plotBarcodeRanks.md)
if drawing that plot is all you want; it calls this reader one run at a
time and never holds more than one raw matrix in memory.

## What is returned, and why it is not Seurat objects

Sparse count matrices, not
[Seurat](https://satijalab.org/seurat/reference/reexports.html) objects.
A raw matrix has one column per barcode in the GEM well – of order 10^5
to 10^6, the overwhelming majority of them empty droplets – and building
a Seurat object means computing and storing per-cell metadata for all of
them. Nothing downstream of a raw matrix wants that; barcode ranking
needs only the column sums. Filtered matrices, which are already
cell-called, are what
[`seuratFrom10x()`](https://mikemartinez99.github.io/rscprep/reference/seuratFrom10x.md)
wraps in Seurat objects.

## `cellranger multi` runs yield one matrix per run, not per sample

This is the one place where the return shape differs from
[`seuratFrom10x()`](https://mikemartinez99.github.io/rscprep/reference/seuratFrom10x.md).
`cellranger multi` writes a single run-level raw matrix at
`outs/multi/count/raw_feature_bc_matrix.h5`, shared by every demuxed
sample in the pool – cell calling happens once for the whole GEM well,
and only afterwards are the called barcodes assigned to samples. There
is therefore no per-sample raw matrix to read, and a `multi = TRUE` row
expands to one entry named after the run directory rather than to N
entries named after the demuxed samples.

Read that way, a barcode rank curve from a pooled run describes the
pool, which is the level at which the knee is a meaningful diagnostic.

## See also

[`plotBarcodeRanks()`](https://mikemartinez99.github.io/rscprep/reference/plotBarcodeRanks.md),
which consumes these matrices one at a time;
[`seuratFrom10x()`](https://mikemartinez99.github.io/rscprep/reference/seuratFrom10x.md)
for the filtered, cell-called matrices.

## Examples

``` r
if (FALSE) { # \dontrun{
# One run
raw <- rawMatrixFrom10x(path = "/path/to/Sample_1/outs", multi = FALSE)
dim(raw[["Sample_1"]])

# A pool -- one matrix for the whole run, not one per demuxed sample
raw <- rawMatrixFrom10x(path = "/path/to/Pool1/outs", multi = TRUE)
names(raw)

# Several runs from a manifest. Watch memory: these are large.
manifest <- data.frame(
  Path  = c("/path/to/Sample_1/outs", "/path/to/Pool1/outs"),
  multi = c(FALSE, TRUE)
)
raw <- rawMatrixFrom10x(manifest = manifest)

# Column sums are the barcode totals a rank plot is built from
totals <- Matrix::colSums(raw[[1]])
} # }
```
