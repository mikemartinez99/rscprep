# Barcode rank plots for unpooled samples

``` r

library(rscprep)
```

A barcode rank plot is the first thing to look at when a sample comes
back with an unexpected cell count. It plots total UMI count against
barcode rank on log-log axes: a sharp cliff between a plateau of real
cells and a long tail of empty droplets means cell calling had an easy
decision to make, and a shallow gradual slope means it did not.

[`plotBarcodeRanks()`](https://mikemartinez99.github.io/rscprep/reference/plotBarcodeRanks.md)
draws that plot. This vignette covers the simplest case – one biological
sample per GEM well, nothing multiplexed, where each curve describes
exactly one sample and needs no caveats.

`DropletUtils` does the ranking and is a `Suggests`, not an `Imports`,
so install it before you start:

``` r

BiocManager::install(c("DropletUtils", "S4Vectors"))
```

## One sample

``` r

ranks <- plotBarcodeRanks(
  path  = "/path/to/Sample_1/outs",
  multi = FALSE
)

ranks$Sample_1$plot
```

`path` is the same `outs/` directory you would hand
[`seuratFrom10x()`](https://mikemartinez99.github.io/rscprep/reference/seuratFrom10x.md),
and the sample is named the same way (after the folder above `outs`).
What differs is which file gets read.

## Why this reads a different file

Cell calling is exactly the step that separates the raw matrix from the
filtered one. The filtered matrix has already had the empty droplets
removed, so the left-hand cliff (the thing the plot exists to show you)
is gone before you start.
[`plotBarcodeRanks()`](https://mikemartinez99.github.io/rscprep/reference/plotBarcodeRanks.md)
therefore reads the **raw** matrix:

|  | [`seuratFrom10x()`](https://mikemartinez99.github.io/rscprep/reference/seuratFrom10x.md) | [`plotBarcodeRanks()`](https://mikemartinez99.github.io/rscprep/reference/plotBarcodeRanks.md) |
|----|----|----|
| Reads | `outs/filtered_feature_bc_matrix.h5` | `outs/raw_feature_bc_matrix.h5` |
| Contains | barcodes called as cells | every barcode in the GEM well |
| Order | 10^3 to 10^4 columns | 10^5 to 10^6 columns |

Nothing is filtered by this function. The knee and the inflection are
drawn so you can compare them against what Cell Ranger actually kept;
neither is a threshold applied to anything.

## “Unpooled” has two shapes

The distinction that matters for reading these curves is whether the GEM
well held one sample or several.

``` r

# A `cellranger count` run
plotBarcodeRanks(path = "/path/to/Sample_1/outs", multi = FALSE)

# A singleplex Flex or 5' run, processed with `cellranger multi`
plotBarcodeRanks(path = "/path/to/Sample_2/outs", multi = TRUE)
```

Both are unpooled, and both give you one curve for one sample. Set
`multi` to match the pipeline that wrote the directory, not to whether
anything was multiplexed: `cellranger multi` writes `per_sample_outs/`
even when it demultiplexed exactly one sample, and the raw matrix moves
with it:

| `multi` | Raw matrix location                         |
|---------|---------------------------------------------|
| `FALSE` | `outs/raw_feature_bc_matrix.h5`             |
| `TRUE`  | `outs/multi/count/raw_feature_bc_matrix.h5` |

If your `cellranger multi` output has been staged or flattened and the
raw matrix sits directly in `outs/`, that is checked as a fallback, so
both layouts resolve either way. `$h5` on the result records which file
was actually read.

The reason this is worth being careful about: for a genuinely pooled
run, `cellranger multi` writes **one** run-level raw matrix shared by
every demuxed sample, because cell calling happens once per GEM well and
only afterwards are barcodes assigned to samples. So a pool gives you
one curve describing the pool, not one curve per sample. When the well
held a single sample, that run-level matrix *is* that sample’s matrix
and the caveat disappears, which is why this vignette can treat both
rows above the same way.

## A batch of samples

Same manifest
[`seuratFrom10x()`](https://mikemartinez99.github.io/rscprep/reference/seuratFrom10x.md)
takes, so one table drives both halves of the workflow:

``` r

manifest <- data.frame(
  Path = c(
    "/path/to/Sample_1/outs",
    "/path/to/Sample_2/outs",
    "/path/to/Sample_3/outs"
  ),
  multi = FALSE
)

ranks <- plotBarcodeRanks(manifest = manifest)

names(ranks)
#> [1] "Sample_1" "Sample_2" "Sample_3"
```

Runs are processed strictly one at a time – read, ranked, reduced to a
curve, released – so peak memory is set by your largest single sample
rather than by the length of the manifest. That matters here more than
it looks: a raw matrix has a column for every barcode in the well, and
holding a dozen of them at once is what makes this kind of loop fall
over.

## What comes back

One element per sample, each a list:

``` r

sample_1 <- ranks$Sample_1

sample_1$plot                 # a ggplot2 object; nothing was printed or saved
head(sample_1$data)           # the plotted curve
sample_1$knee                 # optimistic cutoff, in UMIs
sample_1$inflection           # permissive cutoff, in UMIs
sample_1$n_barcodes           # barcodes in the raw matrix
sample_1$n_above_inflection   # roughly, the cells you would expect
sample_1$h5                   # the file this was computed from
```

`data` is the curve, not the raw matrix: zero-count barcodes are dropped
(they have no place on a log axis, and in a raw matrix they are the
clear majority) and barcodes tying on total are collapsed to one row,
since they share a rank and plot at one point anyway. That turns
millions of rows into a few thousand – small enough to write to CSV,
which the full per-barcode table is not.

## The check this is actually for

`n_above_inflection` is a rough expected cell count derived
independently of Cell Ranger’s cell caller. Comparing it against what
the filtered matrix actually kept is the whole point of drawing these:

``` r

seurat_list <- seuratFrom10x(manifest = manifest)

data.frame(
  expected = vapply(ranks, function(x) x$n_above_inflection, numeric(1)),
  called   = vapply(seurat_list, ncol, integer(1))
)
#>          expected called
#> Sample_1     4914   4821
#> Sample_2     5230   5104
#> Sample_3    11702   6633
```

The first two agree closely, which is what a clean run looks like.
`Sample_3` does not, and that is the signal: go look at its curve. A
shallow slope with no clear cliff usually means ambient RNA, low
viability, or a wetting failure, and the called cell count from that
sample should be treated with suspicion regardless of which side of the
discrepancy you believe.

## Writing figures out

``` r

ranks <- plotBarcodeRanks(manifest = manifest, outdir = "qc/barcode_ranks")
```

Writes `<sample>_barcode_ranks.csv` and `<sample>_barcode_ranks.png` per
sample, creating the directory if needed, and still returns everything
in memory. Sample names come from directory names, so any character that
is not safe in a filename is replaced with `_` in the file stem only –
the list names are left alone. `width`, `height` and `dpi` are passed
through to
[`ggplot2::ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html).

## Tuning the fit

`lower` is the total UMI count at or below which a barcode is assumed
empty. It sets the ambient profile and bounds the spline fit; it does
**not** filter the curve, and barcodes below it are still plotted.

``` r

plotBarcodeRanks(manifest = manifest, lower = 500)$Sample_1$plot
```

Anything else goes to
[`DropletUtils::barcodeRanks()`](https://rdrr.io/pkg/DropletUtils/man/barcodeRanks.html):

``` r

plotBarcodeRanks(
  manifest     = manifest,
  lower        = 100,
  fit.bounds   = NULL,
  exclude.from = 50,
  df           = 20
)
```

Raise `df` for a stiffer spline if the fitted line is chasing noise in
the tail; `exclude.from` drops the highest-ranked barcodes from the fit,
which helps when a handful of enormous barcodes drag the knee upward.

## Changing the figure

The plots are plain ggplots, so nothing is locked in:

``` r

library(ggplot2)

ranks$Sample_1$plot +
  labs(subtitle = "Re-sequenced 2026-08-14") +
  theme_prepro("clean")
```

And `data` is a plain data frame, so a single figure across every sample
is one `rbind` away:

``` r

curves <- do.call(rbind, lapply(names(ranks), function(s) {
  cbind(Sample = s, ranks[[s]]$data)
}))

ggplot(curves, aes(rank, total, colour = Sample)) +
  geom_line() +
  scale_x_log10() +
  scale_y_log10() +
  scale_colour_prepro() +
  labs(x = "Barcode rank", y = "Total UMI count") +
  theme_prepro("airy")
```

Overlaying samples this way is usually more informative than reading the
per-sample figures side by side – a curve that sits below the others, or
turns over earlier, stands out immediately.

## When a sample fails

Two things can go wrong per sample: the raw matrix is missing, or
`barcodeRanks()` cannot fit a curve to it. Both follow the same rule as
[`seuratFrom10x()`](https://mikemartinez99.github.io/rscprep/reference/seuratFrom10x.md).

In manifest mode the sample is skipped with a warning and the rest of
the batch still runs, so the returned list can be shorter than the
manifest:

``` r

ranks <- plotBarcodeRanks(manifest = manifest)
#> ! Raw matrix for "Sample_2" not found at '/path/to/Sample_2/outs/raw_feature_bc_matrix.h5'
#> v Ranked 2 runs of 3.

setdiff(basename(dirname(manifest$Path)), names(ranks))
#> [1] "Sample_2"
```

In `path` mode the same condition is an error instead – there is no
other sample to fall back on, so failing quietly would leave you with an
empty list and no explanation.

If the raw matrix is genuinely absent, the usual cause is that only the
filtered outputs were copied off the sequencing core’s storage.
`raw_feature_bc_matrix.h5` is large and often gets left behind; there is
no way to reconstruct it from what you have, so it has to be fetched.

## Getting the matrices themselves

[`plotBarcodeRanks()`](https://mikemartinez99.github.io/rscprep/reference/plotBarcodeRanks.md)
never hands you a raw matrix – it reads each one, ranks it, and drops
it. If you want the matrices,
[`rawMatrixFrom10x()`](https://mikemartinez99.github.io/rscprep/reference/rawMatrixFrom10x.md)
takes the same arguments and returns them, with no such promise about
memory:

``` r

raw <- rawMatrixFrom10x(path = "/path/to/Sample_1/outs", multi = FALSE)

dim(raw$Sample_1)
#> [1]   36601 6794880

range(Matrix::colSums(raw$Sample_1))
```

It returns count matrices rather than Seurat objects on purpose.
Wrapping a raw matrix in a Seurat object means computing per-cell
metadata for hundreds of thousands of droplets that are about to be
discarded, and barcode ranking only ever needs the column sums.
