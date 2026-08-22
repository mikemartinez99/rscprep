#== Fallback palette: Okabe-Ito ("Wong") colourblind-safe base of 8, extended by
#== walking each hue darker then lighter. Mirrors the construction dittoSeq uses,
#== so plots look the same whether or not dittoSeq is installed.
.prepro_colors <- c(
    "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#666666",
    "#AD7700", "#1C91D4", "#007756", "#D5C711", "#005685", "#A04700", "#B14380", "#4D4D4D",
    "#FFBE2D", "#80C7EF", "#00F6B3", "#F4EB71", "#06A5FF", "#FF8320", "#D99BBD", "#8C8C8C"
)

#' Colourblind-Safe Discrete Palette for Single-Cell QC
#'
#' Returns the package's categorical colour vector, designed for the sample
#' counts single-cell experiments actually produce -- a 10x Flex pool can carry
#' 16 samples, well past the 8 hues most palettes offer.
#'
#' If \pkg{dittoSeq} is installed, [dittoSeq::dittoColors()] supplies the
#' colours: an Okabe-Ito colourblind-safe base of 8, extended to 40 by stepping
#' each hue darker and lighter. Otherwise an equivalent 24-colour vector built
#' the same way is used, so the package has no hard Bioconductor dependency and
#' figures stay reproducible either way.
#'
#' Colours are assigned in fixed order, never recycled by rank -- the same
#' sample keeps the same colour when other samples are filtered out of a plot.
#'
#' @param n Number of colours to return. `NULL` (default) returns the full
#'   vector. If `n` exceeds the number of available colours the palette wraps,
#'   which is worth avoiding: beyond ~24 groups, facet the plot instead.
#'
#' @return A character vector of hex colour codes.
#'
#' @examples
#' prepro_pal(4)
#' length(prepro_pal())
#'
#' @seealso [scale_fill_prepro()] and [scale_colour_prepro()], which apply this
#'   palette to a plot; [theme_prepro()] for the matching theme.
#'
#' @export
prepro_pal <- function(n = NULL) {
    cols <- if (requireNamespace("dittoSeq", quietly = TRUE)) {
        dittoSeq::dittoColors(reps = 1)
    } else {
        .prepro_colors
    }

    if (is.null(n)) {
        return(cols)
    }

    if (!is.numeric(n) || length(n) != 1 || is.na(n) || n < 1) {
        cli::cli_abort("'n' must be a single positive number, or NULL for the full palette.")
    }

    if (n > length(cols)) {
        cli::cli_warn(
            "Requested {n} colours but only {length(cols)} are available; the palette will repeat."
        )
    }

    rep(cols, length.out = n)
}

#' Discrete Colour and Fill Scales for Single-Cell QC Plots
#'
#' Apply the [prepro_pal()] palette to the `fill` or `colour` aesthetic. These
#' are thin wrappers over [ggplot2::scale_fill_manual()] and
#' [ggplot2::scale_colour_manual()], so every argument those accept -- `name`,
#' `labels`, `breaks`, `na.value`, `guide` -- works here too.
#'
#' `scale_color_prepro()` is an alias for `scale_colour_prepro()`.
#'
#' @param ... Passed to [ggplot2::scale_fill_manual()] or
#'   [ggplot2::scale_colour_manual()].
#'
#' @return A \pkg{ggplot2} scale, to be added to a plot with `+`.
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#'
#' ggplot(seurat@meta.data, aes(x = Ident, y = nCount_RNA, fill = Ident)) +
#'     geom_violin() +
#'     scale_y_log10() +
#'     scale_fill_prepro() +
#'     theme_prepro()
#' }
#'
#' @seealso [prepro_pal()] for the underlying colours, [theme_prepro()] for the
#'   matching theme.
#'
#' @importFrom ggplot2 scale_fill_manual scale_colour_manual
#' @export
scale_fill_prepro <- function(...) {
    ggplot2::scale_fill_manual(values = prepro_pal(), ...)
}

#' @rdname scale_fill_prepro
#' @export
scale_colour_prepro <- function(...) {
    ggplot2::scale_colour_manual(values = prepro_pal(), ...)
}

#' @rdname scale_fill_prepro
#' @export
scale_color_prepro <- scale_colour_prepro
