# Publication-Ready ggplot2 Theme for Single-Cell QC Figures

A clean, minimal ggplot2 theme tuned for the QC figures this package
produces – per-sample violins, count distributions, faceted density
grids. Axis titles are always **bold**, secondary text (subtitle,
caption, legend) is set in muted grey, and panel furniture is kept
recessive so the data carries the figure.

## Usage

``` r
theme_prepro(
  style = c("clean", "airy", "frame"),
  base_size = 13,
  base_family = ""
)
```

## Arguments

- style:

  Panel style: one of `"clean"` (default), `"airy"` or `"frame"`. See
  *Details*.

- base_size:

  Base font size in points. All other text is sized relative to this, so
  this is the one knob to turn when a figure is destined for a small
  multi-panel layout. Default `13`.

- base_family:

  Base font family. Defaults to `""` (the device's default sans-serif),
  which keeps figures portable across machines. Set to e.g.
  `"Helvetica"` or `"Arial"` to match a journal's requirements.

## Value

A complete ggplot2 theme object, to be added to a plot with `+`.

## Details

Three panel styles are available. They share identical typography and
differ only in how the panel is framed, so figures built from a mix of
styles still read as one set:

- `"clean"`:

  No gridlines; a hairline L-shaped axis with ticks. The classic journal
  figure. Best for violins and boxplots, where shape matters more than
  reading exact values off the axis.

- `"airy"`:

  Faint horizontal rules only, with no axis lines or ticks – the grid
  stands in for the axis. Best when a reader needs to compare heights
  across samples, such as median UMI count or percent mitochondrial
  reads.

- `"frame"`:

  A boxed panel with a dotted grid and tinted facet strips, and the
  legend moved to the top so wide facet grids keep their horizontal
  room. Best for multi-panel figures where each facet needs a visible
  boundary.

The theme is *complete* in the ggplot2 sense, so it replaces the active
theme rather than adding to it. Override any individual element by
adding a further
[`ggplot2::theme()`](https://ggplot2.tidyverse.org/reference/theme.html)
call after it.

## See also

[`scale_fill_prepro()`](https://mikemartinez99.github.io/rscprep/reference/scale_fill_prepro.md)
and
[`scale_colour_prepro()`](https://mikemartinez99.github.io/rscprep/reference/scale_fill_prepro.md)
for the matching colour scales,
[`prepro_pal()`](https://mikemartinez99.github.io/rscprep/reference/prepro_pal.md)
for the palette itself.

## Examples

``` r
if (FALSE) { # \dontrun{
library(ggplot2)

seurat_list <- seuratFrom10x(path = "/path/to/Pool1/outs", multi = TRUE)
seurat      <- mergeSeurat(seurat_list)

# UMI counts per cell, one violin per sample
ggplot(seurat@meta.data, aes(x = Ident, y = nCount_RNA, fill = Ident)) +
    geom_violin(scale = "width", colour = NA, alpha = 0.85) +
    geom_boxplot(width = 0.11, outlier.shape = NA, fill = "white") +
    scale_y_log10() +
    scale_fill_prepro() +
    labs(title = "UMI counts per cell", x = NULL, y = "nCount_RNA") +
    theme_prepro("clean") +
    theme(legend.position = "none")

# the same data as a faceted density grid
ggplot(seurat@meta.data, aes(x = nCount_RNA, fill = Ident)) +
    geom_density(alpha = 0.55) +
    scale_x_log10() +
    scale_fill_prepro() +
    facet_wrap(~ Ident) +
    theme_prepro("frame")
} # }
```
