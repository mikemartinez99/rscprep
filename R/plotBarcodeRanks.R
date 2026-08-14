#' Barcode Rank ("Knee") Plots from Raw 10x Matrices
#'
#' Draws the standard barcode rank plot for one or more Cell Ranger runs: total
#' UMI count against barcode rank on log-log axes, with the fitted spline and
#' the knee and inflection points [DropletUtils::barcodeRanks()] identifies.
#' Takes the same `manifest` / `path` / `multiplexed` arguments as
#' [seuratFrom10x()], and returns both the plotted curve and the figure for
#' every run.
#'
#' This is the first plot to look at when a sample's cell count comes back
#' unexpected. A sharp cliff between a plateau of high-count barcodes and a
#' long tail of near-empty ones means cell calling had an easy decision to
#' make. A shallow, gradual slope means it did not -- usually ambient RNA, a
#' low-viability sample, or a wetting failure -- and the called cell count from
#' that run should be treated with suspicion.
#'
#' @section Reading the plot:
#'
#' \describe{
#'   \item{Points}{One per distinct total UMI count. Barcodes tying on total
#'     share a rank and plot at one position, so the curve is drawn from a few
#'     thousand points rather than the millions of barcodes behind it.}
#'   \item{Fitted line}{The smooth spline `barcodeRanks()` fits to the curve;
#'     the knee and inflection are derived from its derivatives. Undefined
#'     below `lower`, so it stops short of the right-hand tail.}
#'   \item{Knee}{Where the curve bends most sharply. The optimistic cutoff --
#'     barcodes above it are cells by any reading.}
#'   \item{Inflection}{Where the curve starts its final descent into the
#'     ambient tail. The permissive cutoff, and usually the closer of the two
#'     to what Cell Ranger actually called.}
#' }
#'
#' Neither line is a threshold this function applies. Nothing is filtered here;
#' both are drawn so you can compare them against the cell count
#' [seuratFrom10x()] returns from the same run's filtered matrix.
#'
#' @section Memory:
#'
#' Raw matrices are large, and runs are processed strictly one at a time: each
#' is read, ranked, reduced to its curve, and released before the next is
#' touched. Peak memory is therefore set by the largest single run, not by the
#' size of the manifest. [rawMatrixFrom10x()] is the function to reach for if
#' you genuinely need the matrices themselves, and it makes no such promise.
#'
#' @inheritParams seuratFrom10x
#' @param lower Total UMI count at or below which a barcode is assumed to be
#'   empty. Used by [DropletUtils::barcodeRanks()] to estimate the ambient
#'   profile and to bound the spline fit; it does **not** filter the curve, and
#'   barcodes below it are still plotted. Default `100`.
#' @param outdir Directory to write a CSV of the curve and a PNG of the figure
#'   for each run, named `<run>_barcode_ranks.csv` / `.png`. Created if it does
#'   not exist. `NULL` (default) writes nothing and returns everything in
#'   memory.
#' @param width,height,dpi Passed to [ggplot2::ggsave()] when `outdir` is
#'   supplied. Ignored otherwise.
#' @param ... Further arguments to [DropletUtils::barcodeRanks()], such as
#'   `fit.bounds`, `exclude.from` or `df`.
#'
#' @return A named list with one element per run, each a list of:
#'
#'   \describe{
#'     \item{`data`}{Data frame of the plotted curve: `barcode`, `rank`,
#'       `total` and `fitted`, one row per distinct total, ordered by rank.
#'       Zero-count barcodes are dropped (they have no place on a log axis) and
#'       ties are collapsed to one row, so this is far smaller than the raw
#'       matrix -- small enough to write to CSV, which the full per-barcode
#'       table is not.}
#'     \item{`plot`}{A \pkg{ggplot2} object, using [theme_prepro()]. Nothing is
#'       printed, so it can be assigned, modified with `+`, or saved.}
#'     \item{`knee`, `inflection`}{The two cutoffs, in UMIs.}
#'     \item{`n_barcodes`}{Total barcodes in the raw matrix, before any were
#'       dropped for plotting.}
#'     \item{`n_above_inflection`}{Barcodes at or above the inflection point --
#'       a rough expected cell count to compare against the filtered matrix.}
#'     \item{`h5`}{The raw H5 file the curve was computed from.}
#'   }
#'
#'   Runs whose raw matrix is missing, or whose curve `barcodeRanks()` cannot
#'   fit, are dropped from the manifest-mode result with a warning; in `path`
#'   mode either raises an error instead.
#'
#' @examples
#' \dontrun{
#' # One run, in memory
#' ranks <- plotBarcodeRanks(path = "/path/to/Sample_1/outs", multiplexed = FALSE)
#' ranks$Sample_1$plot
#' head(ranks$Sample_1$data)
#' ranks$Sample_1$knee
#'
#' # Every run in a manifest, with CSVs and PNGs written out
#' manifest <- data.frame(
#'   Path        = c("/path/to/Sample_1/outs", "/path/to/Pool1/outs"),
#'   multiplexed = c(FALSE, TRUE)
#' )
#' ranks <- plotBarcodeRanks(manifest = manifest, outdir = "qc/barcode_ranks")
#'
#' # Compare the expected cell count against what cell calling kept
#' vapply(ranks, function(x) x$n_above_inflection, numeric(1))
#'
#' # The plots are plain ggplots
#' ranks[[1]]$plot + ggplot2::labs(subtitle = "Run 1, re-sequenced")
#'
#' # One faceted figure across every run
#' curves <- do.call(rbind, lapply(names(ranks), function(s) {
#'   cbind(Sample = s, ranks[[s]]$data)
#' }))
#' ggplot2::ggplot(curves, ggplot2::aes(rank, total, colour = Sample)) +
#'   ggplot2::geom_line() +
#'   ggplot2::scale_x_log10() +
#'   ggplot2::scale_y_log10() +
#'   scale_colour_prepro() +
#'   theme_prepro("airy")
#' }
#'
#' @seealso [rawMatrixFrom10x()] for the matrices this reads,
#'   [DropletUtils::barcodeRanks()] for the ranking itself, [theme_prepro()]
#'   for the theme.
#'
#' @importFrom ggplot2 .data
#' @export
plotBarcodeRanks <- function(manifest = NULL, path = NULL, multiplexed = NULL,
                             lower  = 100,
                             outdir = NULL,
                             width  = 6.5,
                             height = 5,
                             dpi    = 300,
                             ...) {

    requireBarcodeRankPkgs()
    validateSourceArgs(manifest, path, multiplexed)
    if (!is.null(manifest)) validateManifest(manifest)

    if (!is.numeric(lower) || length(lower) != 1 || is.na(lower) || lower < 0) {
        cli::cli_abort("'lower' must be a single non-negative number.")
    }
    outdir <- prepareOutdir(outdir)

    #== In path mode there is no other run to fall back on, so anything that
    #== would skip a run aborts instead. Same contract as seuratFrom10x().
    on_problem <- if (is.null(path)) {
        function(msg) cli::cli_warn("{msg}")
    } else {
        function(msg) cli::cli_abort("{msg}", call = NULL)
    }

    runs <- enumerateRawRuns(manifest, path, multiplexed)

    results <- lapply(seq_len(nrow(runs)), function(i) {
        ident <- runs$ident[i]
        h5    <- resolveRawH5(runs$outs_dir[i], runs$multiplexed[i])

        counts <- readOneRaw10x(
            outs_dir    = runs$outs_dir[i],
            ident       = ident,
            multiplexed = runs$multiplexed[i],
            on_missing  = on_problem
        )
        if (is.null(counts)) return(NULL)

        #== Rank while the matrix is in scope, then drop it before moving on to
        #== the next run -- see the Memory section.
        ranked <- tryCatch(
            DropletUtils::barcodeRanks(counts, lower = lower, ...),
            error = function(e) {
                on_problem(cli::format_inline(
                    "{.val {ident}}: barcode ranking failed ({conditionMessage(e)})"
                ))
                NULL
            }
        )
        rm(counts)
        gc(verbose = FALSE, full = FALSE)

        if (is.null(ranked)) return(NULL)

        result <- summariseBarcodeRanks(ranked, ident, h5)
        if (is.null(result)) {
            on_problem(cli::format_inline(
                "{.val {ident}}: no barcode has a non-zero total count"
            ))
            return(NULL)
        }

        result$plot <- barcodeRankPlot(result, ident)
        if (!is.null(outdir)) writeBarcodeRankOutputs(result, ident, outdir,
                                                      width, height, dpi)
        result
    })
    names(results) <- runs$ident

    results <- results[!vapply(results, is.null, logical(1))]

    cli::cli_alert_success(
        "Ranked {length(results)} run{?s} of {nrow(runs)}."
    )

    results
}

#' Check that the Bioconductor packages barcode ranking needs are installed
#'
#' \pkg{DropletUtils} is a Suggests, not an Imports, so that the package
#' installs from CRAN alone; the check is deferred to the one function that
#' actually needs it.
#'
#' @return `invisible(TRUE)`, or aborts.
#'
#' @noRd
requireBarcodeRankPkgs <- function() {
    needed  <- c("DropletUtils", "S4Vectors")
    missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]

    if (length(missing) > 0) {
        install_call <- paste0(
            "BiocManager::install(", deparse(missing, width.cutoff = 500L), ")"
        )
        cli::cli_abort(c(
            "Barcode rank plots require the {.pkg {missing}} package{?s}, which {?is/are} not installed.",
            "i" = "Install with {.code {install_call}}."
        ))
    }
    invisible(TRUE)
}

#' Resolve and create the output directory
#'
#' @param outdir The caller's `outdir`, possibly `NULL`.
#'
#' @return The normalised directory path, or `NULL`.
#'
#' @noRd
prepareOutdir <- function(outdir) {
    if (is.null(outdir)) {
        return(NULL)
    }
    if (!is.character(outdir) || length(outdir) != 1 || is.na(outdir) || outdir == "") {
        cli::cli_abort("'outdir' must be a single directory path, or NULL to write nothing.")
    }

    if (!dir.exists(outdir)) {
        created <- dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
        if (!created) {
            cli::cli_abort("Could not create {.path {outdir}}.")
        }
        cli::cli_alert_info("Created {.path {outdir}}.")
    }

    outdir
}

#' Reduce a barcodeRanks() result to the curve that gets plotted
#'
#' Two reductions, both necessary rather than cosmetic. Zero-total barcodes
#' have no logarithm and cannot appear on a log axis at all -- and in a raw
#' matrix they are the clear majority of columns. Barcodes tying on total share
#' an averaged rank and so plot at one point, meaning all but one of each tie
#' is redundant ink. Together these turn millions of rows into a few thousand,
#' which is the difference between a figure that renders and one that does not.
#'
#' @param ranked The [DropletUtils::barcodeRanks()] result.
#' @param ident Sample name, for messages.
#' @param h5 Path the matrix was read from, recorded in the result.
#'
#' @return A list with `data`, `knee`, `inflection`, `n_barcodes`,
#'   `n_above_inflection` and `h5`, or `NULL` if no barcode has a non-zero
#'   total.
#'
#' @noRd
summariseBarcodeRanks <- function(ranked, ident, h5) {
    curve <- data.frame(
        barcode = if (is.null(rownames(ranked))) NA_character_ else rownames(ranked),
        rank    = as.numeric(ranked$rank),
        total   = as.numeric(ranked$total),
        fitted  = as.numeric(ranked$fitted),
        stringsAsFactors = FALSE
    )

    n_barcodes <- nrow(curve)
    knee       <- S4Vectors::metadata(ranked)$knee
    inflection <- S4Vectors::metadata(ranked)$inflection

    n_above_inflection <- sum(curve$total >= inflection, na.rm = TRUE)

    curve <- curve[!is.na(curve$total) & curve$total > 0, , drop = FALSE]
    if (nrow(curve) == 0) {
        return(NULL)
    }

    curve <- curve[order(curve$rank), , drop = FALSE]
    curve <- curve[!duplicated(curve$rank), , drop = FALSE]
    rownames(curve) <- NULL

    n_empty <- n_barcodes - n_above_inflection
    cli::cli_alert_info(
        "{.val {ident}}: knee at {round(knee)} UMIs, inflection at {round(inflection)} UMIs; {n_above_inflection} barcode{?s} above the inflection ({n_empty} below)."
    )

    list(
        data               = curve,
        knee               = knee,
        inflection         = inflection,
        n_barcodes         = n_barcodes,
        n_above_inflection = n_above_inflection,
        h5                 = h5
    )
}

#== Okabe-Ito hues, matching prepro_pal(). Fixed rather than taken by index
#== from the palette, because these three roles are not interchangeable
#== categories -- knee and inflection must stay visually distinct from each
#== other and from the fitted line no matter how many samples are plotted.
.br_fitted     <- "#D55E00"
.br_knee       <- "#0072B2"
.br_inflection <- "#009E73"

#' Dashed rule plus label for one cutoff
#'
#' The label is anchored to the leftmost data point rather than to the panel
#' edge: `x = -Inf` is the usual way to pin an annotation to the edge, but a
#' log-transformed position scale turns it into `NaN` and drops the annotation
#' without drawing anything.
#'
#' @param value The cutoff, in UMIs. `NA` yields no layers at all.
#' @param label Name of the cutoff, prefixed to the value.
#' @param colour Line and text colour.
#' @param x_at Data-space x position for the label.
#' @param vjust Vertical justification: negative puts the label above the line.
#'
#' @return A list of \pkg{ggplot2} layers, possibly empty.
#'
#' @noRd
cutoffLayers <- function(value, label, colour, x_at, vjust) {
    if (length(value) != 1 || is.na(value) || value <= 0) {
        cli::cli_warn("No {tolower(label)} point could be estimated; the line is omitted.")
        return(list())
    }

    list(
        ggplot2::geom_hline(yintercept = value,
                            colour     = colour,
                            linetype   = "dashed",
                            linewidth  = 0.45),
        ggplot2::annotate("text",
                          x        = x_at,
                          y        = value,
                          label    = paste0(label, ": ", format(round(value), big.mark = ",")),
                          colour   = colour,
                          hjust    = 0,
                          vjust    = vjust,
                          size     = 3.1,
                          fontface = "bold")
    )
}

#' Build the barcode rank figure for one run
#'
#' @param result A `summariseBarcodeRanks()` result.
#' @param ident Sample name, used as the title.
#'
#' @return A \pkg{ggplot2} object.
#'
#' @noRd
barcodeRankPlot <- function(result, ident) {
    curve <- result$data
    point <- pointAesthetics(nrow(curve))

    #== Multiplicative, not additive: x is log scaled, so a fixed offset would
    #== be invisible at rank 10,000 and enormous at rank 1
    label_x <- min(curve$rank) * 1.2

    ggplot2::ggplot(curve, ggplot2::aes(x = .data$rank, y = .data$total)) +
        ggplot2::geom_point(
            colour = "grey45",
            size   = point$size,
            alpha  = point$alpha,
            na.rm  = TRUE
        ) +
        ggplot2::geom_line(
            ggplot2::aes(y = .data$fitted),
            colour    = .br_fitted,
            linewidth = 0.7,
            na.rm     = TRUE
        ) +
        #== Knee is always the higher of the two cutoffs, so its label goes
        #== above its line and the inflection's below its own; the pair stays
        #== legible even when the two nearly coincide.
        cutoffLayers(result$knee, "Knee", .br_knee, label_x, vjust = -0.55) +
        cutoffLayers(result$inflection, "Inflection", .br_inflection,
                     label_x, vjust = 1.35) +
        ggplot2::scale_x_log10() +
        ggplot2::scale_y_log10() +
        ggplot2::labs(
            title    = ident,
            subtitle = paste0(
                format(result$n_above_inflection, big.mark = ","),
                " barcodes above the inflection point of ",
                format(result$n_barcodes, big.mark = ","), " raw"
            ),
            x = "Barcode rank",
            y = "Total UMI count"
        ) +
        theme_prepro("airy")
}

#' Write one run's curve and figure to disk
#'
#' @param result A `summariseBarcodeRanks()` result, with `plot` attached.
#' @param ident Sample name, used for the filenames.
#' @param outdir Directory to write into.
#' @param width,height,dpi Passed to [ggplot2::ggsave()].
#'
#' @return `invisible(NULL)`.
#'
#' @noRd
writeBarcodeRankOutputs <- function(result, ident, outdir, width, height, dpi) {
    #== Sample names come from directory names, which are not guaranteed to be
    #== safe filename stems on every platform
    stem <- gsub("[^A-Za-z0-9._-]+", "_", ident)

    csv_path <- file.path(outdir, paste0(stem, "_barcode_ranks.csv"))
    png_path <- file.path(outdir, paste0(stem, "_barcode_ranks.png"))

    utils::write.csv(result$data, file = csv_path, row.names = FALSE)
    ggplot2::ggsave(png_path, plot = result$plot,
                    width = width, height = height, dpi = dpi)

    cli::cli_alert_success("{.val {ident}}: wrote {.path {csv_path}} and {.path {png_path}}.")
    invisible(NULL)
}
