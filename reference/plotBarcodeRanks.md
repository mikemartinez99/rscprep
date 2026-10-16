# Barcode Rank ("Knee") Plots from Raw 10x Matrices

Draws the standard barcode rank plot for one or more Cell Ranger runs:
total UMI count against barcode rank on log-log axes, with the fitted
spline and the knee and inflection points
[`DropletUtils::barcodeRanks()`](https://rdrr.io/pkg/DropletUtils/man/barcodeRanks.html)
identifies. Takes the same `manifest` / `path` / `multi` arguments as
[`seuratFrom10x()`](https://mikemartinez99.github.io/rscprep/reference/seuratFrom10x.md),
and returns both the plotted curve and the figure for every run.

## Usage

``` r
plotBarcodeRanks(
  manifest = NULL,
  path = NULL,
  multi = NULL,
  lower = 100,
  outdir = NULL,
  width = 6.5,
  height = 5,
  dpi = 300,
  ...,
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

- lower:

  Total UMI count at or below which a barcode is assumed to be empty.
  Used by
  [`DropletUtils::barcodeRanks()`](https://rdrr.io/pkg/DropletUtils/man/barcodeRanks.html)
  to estimate the ambient profile and to bound the spline fit; it does
  **not** filter the curve, and barcodes below it are still plotted.
  Default `100`.

- outdir:

  Directory to write a CSV of the curve and a PNG of the figure for each
  run, named `<run>_barcode_ranks.csv` / `.png`. Created if it does not
  exist. `NULL` (default) writes nothing and returns everything in
  memory.

- width, height, dpi:

  Passed to
  [`ggplot2::ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html)
  when `outdir` is supplied. Ignored otherwise.

- ...:

  Further arguments to
  [`DropletUtils::barcodeRanks()`](https://rdrr.io/pkg/DropletUtils/man/barcodeRanks.html),
  such as `fit.bounds`, `exclude.from` or `df`.

- multiplexed:

  Deprecated. The former name of `multi`; still accepted, with a
  warning.

## Value

A named list with one element per run, each a list of:

- `data`:

  Data frame of the plotted curve: `barcode`, `rank`, `total` and
  `fitted`, one row per distinct total, ordered by rank. Zero-count
  barcodes are dropped (they have no place on a log axis) and ties are
  collapsed to one row, so this is far smaller than the raw matrix –
  small enough to write to CSV, which the full per-barcode table is not.

- `plot`:

  A ggplot2 object, using
  [`theme_prepro()`](https://mikemartinez99.github.io/rscprep/reference/theme_prepro.md).
  Nothing is printed, so it can be assigned, modified with `+`, or
  saved.

- `knee`, `inflection`:

  The two cutoffs, in UMIs.

- `n_barcodes`:

  Total barcodes in the raw matrix, before any were dropped for
  plotting.

- `n_above_inflection`:

  Barcodes at or above the inflection point – a rough expected cell
  count to compare against the filtered matrix.

- `h5`:

  The raw H5 file the curve was computed from.

Runs whose raw matrix is missing, or whose curve `barcodeRanks()` cannot
fit, are dropped from the manifest-mode result with a warning; in `path`
mode either raises an error instead.

## Details

This is the first plot to look at when a sample's cell count comes back
unexpected. A sharp cliff between a plateau of high-count barcodes and a
long tail of near-empty ones means cell calling had an easy decision to
make. A shallow, gradual slope means it did not – usually ambient RNA, a
low-viability sample, or a wetting failure – and the called cell count
from that run should be treated with suspicion.

## Reading the plot

- Points:

  One per distinct total UMI count. Barcodes tying on total share a rank
  and plot at one position, so the curve is drawn from a few thousand
  points rather than the millions of barcodes behind it.

- Fitted line:

  The smooth spline `barcodeRanks()` fits to the curve; the knee and
  inflection are derived from its derivatives. Undefined below `lower`,
  so it stops short of the right-hand tail.

- Knee:

  Where the curve bends most sharply. The optimistic cutoff – barcodes
  above it are cells by any reading.

- Inflection:

  Where the curve starts its final descent into the ambient tail. The
  permissive cutoff, and usually the closer of the two to what Cell
  Ranger actually called.

Neither line is a threshold this function applies. Nothing is filtered
here; both are drawn so you can compare them against the cell count
[`seuratFrom10x()`](https://mikemartinez99.github.io/rscprep/reference/seuratFrom10x.md)
returns from the same run's filtered matrix.

## Memory

Raw matrices are large, and runs are processed strictly one at a time:
each is read, ranked, reduced to its curve, and released before the next
is touched. Peak memory is therefore set by the largest single run, not
by the size of the manifest.
[`rawMatrixFrom10x()`](https://mikemartinez99.github.io/rscprep/reference/rawMatrixFrom10x.md)
is the function to reach for if you genuinely need the matrices
themselves, and it makes no such promise.

## See also

[`rawMatrixFrom10x()`](https://mikemartinez99.github.io/rscprep/reference/rawMatrixFrom10x.md)
for the matrices this reads,
[`DropletUtils::barcodeRanks()`](https://rdrr.io/pkg/DropletUtils/man/barcodeRanks.html)
for the ranking itself,
[`theme_prepro()`](https://mikemartinez99.github.io/rscprep/reference/theme_prepro.md)
for the theme.

## Examples

``` r
if (FALSE) { # \dontrun{
# One run, in memory
ranks <- plotBarcodeRanks(path = "/path/to/Sample_1/outs", multi = FALSE)
ranks$Sample_1$plot
head(ranks$Sample_1$data)
ranks$Sample_1$knee

# Every run in a manifest, with CSVs and PNGs written out
manifest <- data.frame(
  Path  = c("/path/to/Sample_1/outs", "/path/to/Pool1/outs"),
  multi = c(FALSE, TRUE)
)
ranks <- plotBarcodeRanks(manifest = manifest, outdir = "qc/barcode_ranks")

# Compare the expected cell count against what cell calling kept
vapply(ranks, function(x) x$n_above_inflection, numeric(1))

# The plots are plain ggplots
ranks[[1]]$plot + ggplot2::labs(subtitle = "Run 1, re-sequenced")

# One faceted figure across every run
curves <- do.call(rbind, lapply(names(ranks), function(s) {
  cbind(Sample = s, ranks[[s]]$data)
}))
ggplot2::ggplot(curves, ggplot2::aes(rank, total, colour = Sample)) +
  ggplot2::geom_line() +
  ggplot2::scale_x_log10() +
  ggplot2::scale_y_log10() +
  scale_colour_prepro() +
  theme_prepro("airy")
} # }
```
