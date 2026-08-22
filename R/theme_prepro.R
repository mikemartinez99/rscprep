#' Publication-Ready ggplot2 Theme for Single-Cell QC Figures
#'
#' A clean, minimal \pkg{ggplot2} theme tuned for the QC figures this package
#' produces -- per-sample violins, count distributions, faceted density grids.
#' Axis titles are always **bold**, secondary text (subtitle, caption, legend)
#' is set in muted grey, and panel furniture is kept recessive so the data
#' carries the figure.
#'
#' Three panel styles are available. They share identical typography and differ
#' only in how the panel is framed, so figures built from a mix of styles still
#' read as one set:
#'
#' \describe{
#'   \item{`"clean"`}{No gridlines; a hairline L-shaped axis with ticks. The
#'     classic journal figure. Best for violins and boxplots, where shape
#'     matters more than reading exact values off the axis.}
#'   \item{`"airy"`}{Faint horizontal rules only, with no axis lines or ticks --
#'     the grid stands in for the axis. Best when a reader needs to compare
#'     heights across samples, such as median UMI count or percent
#'     mitochondrial reads.}
#'   \item{`"frame"`}{A boxed panel with a dotted grid and tinted facet strips,
#'     and the legend moved to the top so wide facet grids keep their
#'     horizontal room. Best for multi-panel figures where each facet needs a
#'     visible boundary.}
#' }
#'
#' The theme is *complete* in the \pkg{ggplot2} sense, so it replaces the
#' active theme rather than adding to it. Override any individual element by
#' adding a further [ggplot2::theme()] call after it.
#'
#' @param style Panel style: one of `"clean"` (default), `"airy"` or
#'   `"frame"`. See *Details*.
#' @param base_size Base font size in points. All other text is sized relative
#'   to this, so this is the one knob to turn when a figure is destined for a
#'   small multi-panel layout. Default `13`.
#' @param base_family Base font family. Defaults to `""` (the device's default
#'   sans-serif), which keeps figures portable across machines. Set to e.g.
#'   `"Helvetica"` or `"Arial"` to match a journal's requirements.
#'
#' @return A complete \pkg{ggplot2} theme object, to be added to a plot with
#'   `+`.
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#'
#' seurat_list <- seuratFrom10x(path = "/path/to/Pool1/outs", multi = TRUE)
#' seurat      <- mergeSeurat(seurat_list)
#'
#' # UMI counts per cell, one violin per sample
#' ggplot(seurat@meta.data, aes(x = Ident, y = nCount_RNA, fill = Ident)) +
#'     geom_violin(scale = "width", colour = NA, alpha = 0.85) +
#'     geom_boxplot(width = 0.11, outlier.shape = NA, fill = "white") +
#'     scale_y_log10() +
#'     scale_fill_prepro() +
#'     labs(title = "UMI counts per cell", x = NULL, y = "nCount_RNA") +
#'     theme_prepro("clean") +
#'     theme(legend.position = "none")
#'
#' # the same data as a faceted density grid
#' ggplot(seurat@meta.data, aes(x = nCount_RNA, fill = Ident)) +
#'     geom_density(alpha = 0.55) +
#'     scale_x_log10() +
#'     scale_fill_prepro() +
#'     facet_wrap(~ Ident) +
#'     theme_prepro("frame")
#' }
#'
#' @seealso [scale_fill_prepro()] and [scale_colour_prepro()] for the matching
#'   colour scales, [prepro_pal()] for the palette itself.
#'
#' @importFrom ggplot2 theme element_text element_line element_rect
#'   element_blank margin rel
#' @export
theme_prepro <- function(style = c("clean", "airy", "frame"),
                         base_size   = 13,
                         base_family = "") {
    style <- match.arg(style)

    if (!is.numeric(base_size) || length(base_size) != 1 || is.na(base_size) || base_size <= 0) {
        cli::cli_abort("'base_size' must be a single positive number.")
    }

    base <- preproThemeBase(base_size = base_size, base_family = base_family)

    #== Only the panel furniture differs between styles; typography is shared
    panel <- switch(
        style,
        clean = ggplot2::theme(
            panel.grid        = ggplot2::element_blank(),
            axis.line         = ggplot2::element_line(colour = "grey20", linewidth = 0.4),
            axis.ticks        = ggplot2::element_line(colour = "grey20", linewidth = 0.4),
            axis.ticks.length = grid::unit(3.5, "pt"),
            strip.background  = ggplot2::element_blank()
        ),
        airy = ggplot2::theme(
            panel.grid.major.y = ggplot2::element_line(colour = "grey88", linewidth = 0.35),
            panel.grid.major.x = ggplot2::element_blank(),
            panel.grid.minor   = ggplot2::element_blank(),
            axis.line          = ggplot2::element_blank(),
            axis.ticks         = ggplot2::element_blank(),
            axis.text          = ggplot2::element_text(colour = "grey35"),
            strip.background   = ggplot2::element_blank()
        ),
        frame = ggplot2::theme(
            panel.border       = ggplot2::element_rect(colour = "grey25", fill = NA, linewidth = 0.5),
            panel.grid.major   = ggplot2::element_line(colour = "grey90", linewidth = 0.3,
                                                       linetype = "dotted"),
            panel.grid.minor   = ggplot2::element_blank(),
            axis.ticks         = ggplot2::element_line(colour = "grey25", linewidth = 0.4),
            axis.ticks.length  = grid::unit(3, "pt"),
            strip.background   = ggplot2::element_rect(fill = "grey95", colour = NA),
            panel.spacing      = grid::unit(0.9, "lines"),
            legend.position    = "top",
            legend.justification = "left"
        )
    )

    base + panel
}

#' Shared typographic base for [theme_prepro()]
#'
#' Everything the three panel styles agree on: bold axis titles, muted grey
#' secondary text, and title/caption blocks aligned to the plot edge rather
#' than the panel. Split out so the styles can differ in panel furniture only.
#'
#' @inheritParams theme_prepro
#'
#' @return An incomplete \pkg{ggplot2} theme, built on
#'   [ggplot2::theme_minimal()].
#'
#' @importFrom ggplot2 theme_minimal
#' @noRd
preproThemeBase <- function(base_size = 13, base_family = "") {
    ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
        #== The ask: axis titles always bold
        axis.title            = ggplot2::element_text(face = "bold", colour = "grey15"),
        axis.title.x          = ggplot2::element_text(margin = ggplot2::margin(t = 9)),
        axis.title.y          = ggplot2::element_text(margin = ggplot2::margin(r = 9)),
        axis.text             = ggplot2::element_text(colour = "grey30"),

        #== Title block aligned to the plot edge, not the y-axis
        plot.title            = ggplot2::element_text(face = "bold", size = ggplot2::rel(1.15),
                                                      colour = "grey10",
                                                      margin = ggplot2::margin(b = 3)),
        plot.subtitle         = ggplot2::element_text(size = ggplot2::rel(0.92), colour = "grey40",
                                                      margin = ggplot2::margin(b = 11)),
        plot.caption          = ggplot2::element_text(size = ggplot2::rel(0.78), colour = "grey55",
                                                      hjust = 0,
                                                      margin = ggplot2::margin(t = 11)),
        plot.title.position   = "plot",
        plot.caption.position = "plot",

        legend.title          = ggplot2::element_text(face = "bold", size = ggplot2::rel(0.88)),
        legend.text           = ggplot2::element_text(size = ggplot2::rel(0.82), colour = "grey25"),
        legend.key.size       = grid::unit(0.85, "lines"),
        legend.background     = ggplot2::element_blank(),

        strip.text            = ggplot2::element_text(face = "bold", size = ggplot2::rel(0.92),
                                                      colour = "grey15",
                                                      margin = ggplot2::margin(5, 5, 5, 5)),
        plot.margin           = ggplot2::margin(13, 15, 11, 13)
    )
}
