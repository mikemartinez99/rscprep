# Discrete Colour and Fill Scales for Single-Cell QC Plots

Apply the
[`prepro_pal()`](https://mikemartinez99.github.io/rscprep/reference/prepro_pal.md)
palette to the `fill` or `colour` aesthetic. These are thin wrappers
over
[`ggplot2::scale_fill_manual()`](https://ggplot2.tidyverse.org/reference/scale_manual.html)
and
[`ggplot2::scale_colour_manual()`](https://ggplot2.tidyverse.org/reference/scale_manual.html),
so every argument those accept – `name`, `labels`, `breaks`, `na.value`,
`guide` – works here too.

## Usage

``` r
scale_fill_prepro(...)

scale_colour_prepro(...)

scale_color_prepro(...)
```

## Arguments

- ...:

  Passed to
  [`ggplot2::scale_fill_manual()`](https://ggplot2.tidyverse.org/reference/scale_manual.html)
  or
  [`ggplot2::scale_colour_manual()`](https://ggplot2.tidyverse.org/reference/scale_manual.html).

## Value

A ggplot2 scale, to be added to a plot with `+`.

## Details

`scale_color_prepro()` is an alias for `scale_colour_prepro()`.

## See also

[`prepro_pal()`](https://mikemartinez99.github.io/rscprep/reference/prepro_pal.md)
for the underlying colours,
[`theme_prepro()`](https://mikemartinez99.github.io/rscprep/reference/theme_prepro.md)
for the matching theme.

## Examples

``` r
if (FALSE) { # \dontrun{
library(ggplot2)

ggplot(seurat@meta.data, aes(x = Ident, y = nCount_RNA, fill = Ident)) +
    geom_violin() +
    scale_y_log10() +
    scale_fill_prepro() +
    theme_prepro()
} # }
```
