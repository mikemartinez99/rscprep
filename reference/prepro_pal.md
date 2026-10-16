# Colourblind-Safe Discrete Palette for Single-Cell QC

Returns the package's categorical colour vector, designed for the sample
counts single-cell experiments actually produce – a 10x Flex pool can
carry 16 samples, well past the 8 hues most palettes offer.

## Usage

``` r
prepro_pal(n = NULL)
```

## Arguments

- n:

  Number of colours to return. `NULL` (default) returns the full vector.
  If `n` exceeds the number of available colours the palette wraps,
  which is worth avoiding: beyond ~24 groups, facet the plot instead.

## Value

A character vector of hex colour codes.

## Details

If dittoSeq is installed,
[`dittoSeq::dittoColors()`](https://rdrr.io/pkg/dittoSeq/man/dittoColors.html)
supplies the colours: an Okabe-Ito colourblind-safe base of 8, extended
to 40 by stepping each hue darker and lighter. Otherwise an equivalent
24-colour vector built the same way is used, so the package has no hard
Bioconductor dependency and figures stay reproducible either way.

Colours are assigned in fixed order, never recycled by rank – the same
sample keeps the same colour when other samples are filtered out of a
plot.

## See also

[`scale_fill_prepro()`](https://mikemartinez99.github.io/rscprep/reference/scale_fill_prepro.md)
and
[`scale_colour_prepro()`](https://mikemartinez99.github.io/rscprep/reference/scale_fill_prepro.md),
which apply this palette to a plot;
[`theme_prepro()`](https://mikemartinez99.github.io/rscprep/reference/theme_prepro.md)
for the matching theme.

## Examples

``` r
prepro_pal(4)
#> [1] "#E69F00" "#56B4E9" "#009E73" "#F0E442"
length(prepro_pal())
#> [1] 40
```
