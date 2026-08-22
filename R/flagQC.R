#' Flag Outlier Cells by Median Absolute Deviation
#'
#' Marks cells whose QC metrics sit an unusual distance from the rest of the
#' data, using the median absolute deviation (MAD) rather than a fixed cutoff.
#' For each metric a centre and spread are estimated, a cell is flagged if it
#' falls more than `nmads` MADs from the centre, and the result is written back
#' as a logical metadata column. Nothing is subset or removed -- flagging and
#' filtering are kept separate so the flags can be inspected, and
#' cross-tabulated against sample or cell type, before any cells are dropped.
#'
#' MAD-based thresholds adapt to the data. A hard cutoff carried between
#' experiments discards good cells from a low-depth run and keeps bad ones from
#' a deep run; a threshold derived from the data's own distribution does
#' neither. The cost is that it assumes most cells it is estimated from are
#' good -- the median and MAD are computed over all of them -- so on a badly
#' failed sample it will happily call the failure normal and flag the
#' survivors.
#'
#' @section Where the thresholds come from:
#'
#' `group_by` decides what "the rest of the data" means, and it is the most
#' consequential argument here:
#'
#' \describe{
#'   \item{`NULL` (default)}{One global threshold per metric, estimated across
#'     every cell in the object. Right for a single sample, and for a merged
#'     object whose samples were captured and sequenced together.}
#'   \item{`"Ident"`}{A separate threshold per sample --
#'     [seuratFrom10x()] sets `Ident` to the sample name. Right whenever
#'     sequencing depth or cell recovery varies between samples, which is the
#'     usual case. A global cutoff on such an object flags most of the
#'     shallowest sample and almost none of the deepest, regardless of whether
#'     either actually failed.}
#'   \item{any other metadata column}{A separate threshold per level of that
#'     column: a capture batch, a 10x lane, a tissue. Group by the technical
#'     variable you want the threshold to absorb.}
#' }
#'
#' Grouping is not free. Each group is judged only against itself, so a
#' uniformly bad sample looks fine from the inside and a group with few cells
#' gives an unstable estimate. Groups with fewer than three cells carrying the
#' metric, or where over half the cells share one value, are skipped rather
#' than given a meaningless threshold; they are named in a warning and appear
#' in the summary table with `NA` thresholds. The two views are worth
#' comparing -- if per-sample and global flags disagree sharply, the
#' disagreement is itself the finding.
#'
#' @section Choosing a bound:
#'
#' Which tail matters depends on the metric, so the bound is set per metric
#' rather than globally:
#'
#' \describe{
#'   \item{`"lower"`}{Only cells below the threshold. The right choice for
#'     `nFeature_RNA` and `nCount_RNA`, where the low tail is empty droplets
#'     and dying cells but the high tail is largely real, well-captured cells.}
#'   \item{`"higher"`}{Only cells above it. The right choice for `pctMito`,
#'     where high is bad and low is not.}
#'   \item{`"both"`}{Either tail. Useful for count metrics when doublets are
#'     the concern, since the high tail is where they sit.}
#' }
#'
#' @section The log transform:
#'
#' Count metrics are right-skewed, which drags the median up and inflates the
#' MAD, so `log = TRUE` (the default) estimates the centre and spread on a
#' `log2` scale. Thresholds are transformed back before they are applied,
#' reported and stored, so the values in the summary table are always on the
#' metric's own scale and can be compared against it directly. If a metric
#' contains zeroes the transform becomes `log2(x + 1)` and the back-transform
#' compensates; negative values cannot be log-transformed at all, and are an
#' error rather than a silent `NaN`.
#'
#' Set `log = FALSE` for metrics that are already roughly symmetric, such as
#' `log10GenesPerUMI` or a percentage that does not pile up against zero.
#'
#' @param obj A single [Seurat][Seurat::CreateSeuratObject] object carrying the
#'   metrics to be flagged, typically the output of [diagnoseQC()].
#' @param metrics Metadata columns to flag on. Either an unnamed character
#'   vector of column names, which all take the function-level `bound`, or a
#'   named vector of `column = bound` pairs, which set the bound per metric.
#'   The two forms can be mixed: `c("nCount_RNA", pctMito = "higher")` flags
#'   counts using `bound` and mitochondrial content on the high side only.
#' @param nmads Number of MADs from the median at which a cell becomes an
#'   outlier. Must be a single positive number. `3` (the default) is the
#'   conventional choice; `2.5` is more aggressive, `5` keeps all but the
#'   extremes.
#' @param bound Default tail to flag, applied to every metric not given its own
#'   bound in `metrics`. One of `"both"` (the default), `"lower"` or
#'   `"higher"`. See *Choosing a bound*.
#' @param log Estimate the median and MAD on a `log2` scale. `TRUE` by default,
#'   which is what count metrics want. See *The log transform*.
#' @param group_by Discrete metadata column within which to estimate
#'   thresholds, so that each cell is judged against its own group rather than
#'   against the whole object. `NULL` (default) estimates one global threshold
#'   per metric; `"Ident"` estimates one per sample. See *Where the thresholds
#'   come from*.
#'
#' @return The input object, with one logical metadata column added per metric,
#'   named `<metric>_<nmads>_<bound>_outliers`. `TRUE` marks an outlier;
#'   cells whose metric is `NA`, whose `group_by` value is `NA`, or whose group
#'   was skipped stay `NA`. The column name does not record `group_by`, so
#'   re-running the same metric with different grouping overwrites the previous
#'   flags rather than sitting alongside them.
#'
#'   A summary is stored in `obj@misc$flagQC` as a data frame, one row per
#'   metric per group, with columns:
#'
#'   \describe{
#'     \item{`metric`, `bound`, `nmads`, `scale`}{What was asked for, and
#'       whether the estimate was made on the `"log2"` or `"identity"` scale.}
#'     \item{`group_by`, `group`}{The column grouped on and the level this row
#'       describes, both `NA` when thresholds were estimated globally.}
#'     \item{`center`, `lower`, `upper`}{The median and the two thresholds, on
#'       the metric's own scale. The threshold for a tail that was not flagged
#'       is `NA`, not an unenforced number.}
#'     \item{`n_lower`, `n_upper`}{Cells flagged by each threshold, again `NA`
#'       for a tail that was not flagged. The two are disjoint, so where both
#'       apply they sum to `n_flagged`.}
#'     \item{`n_flagged`, `pct_flagged`, `n_missing`, `n_cells`}{Totals within
#'       the group, where `n_missing` counts cells the metric was `NA` for and
#'       so could not be judged either way.}
#'     \item{`column`}{The metadata column written, or `NA` for a group that
#'       was skipped.}
#'   }
#'
#'   Re-running on the same object appends to that table rather than replacing
#'   it, so it stays a record of every flag column present. Rows describing a
#'   column that is being overwritten are replaced.
#'
#' @examples
#' \dontrun{
#' seurat <- diagnoseQC(seurat, organism = "mmu")
#'
#' # one global threshold per metric
#' seurat <- flagQC(seurat, c(nCount_RNA = "lower", pctMito = "higher"))
#'
#' # the same, but each sample judged against itself
#' seurat <- flagQC(
#'     seurat,
#'     metrics  = c(nCount_RNA = "lower", nFeature_RNA = "lower", pctMito = "higher"),
#'     group_by = "Ident"
#' )
#'
#' # per-sample thresholds and counts
#' seurat@misc$flagQC[, c("metric", "group", "lower", "upper", "n_flagged")]
#'
#' # which samples lose the most cells?
#' qc <- seurat@misc$flagQC
#' with(qc, tapply(pct_flagged, group, max))
#'
#' # group by a technical variable instead
#' seurat <- flagQC(seurat, "nCount_RNA", bound = "lower", group_by = "capture_batch")
#'
#' # nothing has been removed yet -- look first
#' table(seurat$Ident, seurat$pctMito_3_higher_outliers)
#'
#' # then filter on whichever flags you trust
#' flags <- grep("_outliers$", colnames(seurat@meta.data), value = TRUE)
#' keep  <- !Reduce(`|`, seurat@meta.data[flags])
#' seurat <- seurat[, which(keep)]
#' }
#'
#' @seealso [diagnoseQC()], which calculates the metrics this flags on.
#'
#' @export
flagQC <- function(obj,
                   metrics  = c(nFeature_RNA = "both", nCount_RNA = "both"),
                   nmads    = 3,
                   bound    = c("both", "lower", "higher"),
                   log      = TRUE,
                   group_by = NULL) {

    arg_name <- deparse(substitute(obj))
    validateQCObject(obj, arg_name, fn = "flagQC")

    bound  <- match.arg(bound)
    bounds <- c("lower", "higher", "both")
    md     <- obj@meta.data

    if (!is.numeric(nmads) || length(nmads) != 1 || is.na(nmads) || nmads <= 0) {
        cli::cli_abort("{.arg nmads} must be a single positive number.")
    }
    if (!is.logical(log) || length(log) != 1 || is.na(log)) {
        cli::cli_abort("{.arg log} must be either {.val {TRUE}} or {.val {FALSE}}.")
    }
    if (length(metrics) == 0) {
        cli::cli_abort("{.arg metrics} is empty -- provide at least one column.")
    }
    if (!is.character(metrics)) {
        cli::cli_abort(
            "{.arg metrics} must be a character vector, not {.cls {class(metrics)[1]}}."
        )
    }

    #== Unnamed entries name a column and inherit the function-level `bound`;
    #== named entries are `column = bound` and override it.
    nms <- names(metrics)
    if (is.null(nms)) {
        metric_names  <- unname(metrics)
        metric_bounds <- rep(bound, length(metrics))
    } else {
        unnamed       <- nms == "" | is.na(nms)
        metric_names  <- ifelse(unnamed, unname(metrics), nms)
        metric_bounds <- ifelse(unnamed, bound, unname(metrics))
    }

    bad <- setdiff(metric_bounds, bounds)
    if (length(bad) > 0) {
        cli::cli_abort(c(
            "Invalid bound{?s} in {.arg metrics}: {.val {bad}}.",
            "i" = "Each bound must be one of {.val {bounds}}.",
            "i" = "Entries are {.code column = bound}, so a bare
                   {.code c(\"nCount_RNA\")} takes the function-level
                   {.arg bound} instead."
        ))
    }

    absent <- setdiff(metric_names, colnames(md))
    if (length(absent) > 0) {
        cli::cli_abort(c(
            "{.arg metrics} names {length(absent)} column{?s} not in
             {.arg {arg_name}}: {.val {absent}}.",
            "i" = "Percentage and complexity metrics are added by
                   {.fn diagnoseQC}."
        ))
    }

    #== Two entries for the same metric and bound would write the same column
    #== twice, with only the second surviving.
    combos <- paste(metric_names, metric_bounds)
    if (anyDuplicated(combos)) {
        dupes <- unique(metric_names[duplicated(combos)])
        cli::cli_abort(
            "{.arg metrics} repeats the same metric and bound: {.val {dupes}}."
        )
    }

    n_cells <- ncol(obj)
    groups  <- resolveFlagGroups(group_by, md, arg_name, n_cells)
    levs    <- levels(groups)

    cli::cli_alert_info(
        "Flagging outliers in {length(metric_names)} metric{?s} at {nmads} MAD{?s} \\
         from the median{if (log) ' (log2 scale)' else ''}, \\
         {if (is.null(group_by)) 'across all cells' else \
          cli::format_inline('within {length(levs)} group{?s} of {.field {group_by}}')}."
    )

    rows <- list()

    for (k in seq_along(metric_names)) {
        i <- metric_names[k]
        b <- metric_bounds[k]
        x <- md[[i]]

        if (!is.numeric(x)) {
            cli::cli_abort(
                "Metric {.field {i}} is not numeric, but {.cls {class(x)[1]}}."
            )
        }
        if (log && any(x < 0, na.rm = TRUE)) {
            cli::cli_abort(c(
                "Metric {.field {i}} has negative values and cannot be log-transformed.",
                "i" = "Pass {.code log = FALSE} to estimate on the raw scale."
            ))
        }

        use_lo <- b %in% c("lower", "both")
        use_hi <- b %in% c("higher", "both")
        column <- paste0(i, "_", nmads, "_", b, "_outliers")

        #== NA until a cell is judged: cells with no metric, no group, or in a
        #== group too degenerate to estimate from are never set to FALSE, so an
        #== unjudged cell is never mistaken for one that passed.
        outliers    <- rep(NA, n_cells)
        metric_rows <- vector("list", length(levs))
        skipped     <- character(0)
        reasons     <- character(0)

        for (g in seq_along(levs)) {
            idx <- which(!is.na(groups) & groups == levs[g])
            xg  <- x[idx]
            thr <- madThresholds(xg, nmads, log)

            n_missing <- sum(is.na(xg))

            if (!thr$ok) {
                skipped <- c(skipped, levs[g])
                reasons <- c(reasons, thr$reason)
                metric_rows[[g]] <- flagQCRow(
                    i, b, nmads, log, group_by, levs[g],
                    NA_real_, NA_real_, NA_real_,
                    NA_integer_, NA_integer_, NA_integer_,
                    n_missing, length(idx), NA_character_
                )
                next
            }

            #== log2() is monotonic, so comparing the raw metric against
            #== back-transformed thresholds flags exactly the cells the
            #== reported thresholds say it should.
            below <- !is.na(xg) & xg < thr$lower
            above <- !is.na(xg) & xg > thr$upper

            judged <- rep(FALSE, length(idx))
            if (use_lo) judged <- judged | below
            if (use_hi) judged <- judged | above
            judged[is.na(xg)] <- NA
            outliers[idx] <- judged

            metric_rows[[g]] <- flagQCRow(
                i, b, nmads, log, group_by, levs[g],
                thr$center,
                if (use_lo) thr$lower else NA_real_,
                if (use_hi) thr$upper else NA_real_,
                if (use_lo) sum(below) else NA_integer_,
                if (use_hi) sum(above) else NA_integer_,
                sum(judged, na.rm = TRUE),
                n_missing, length(idx), column
            )
        }

        metric_df <- do.call(rbind, metric_rows)

        #== One warning per metric naming every skipped group, rather than one
        #== warning per group -- a 20-sample object would bury the console.
        if (length(skipped) > 0 && is.null(group_by)) {
            cli::cli_warn(c(
                "!" = "{.field {i}}: no threshold could be estimated ({reasons}).",
                "x" = "Every cell's flag is left as {.val {NA}}, not {.val {FALSE}}."
            ))
        } else if (length(skipped) > 0) {
            detail <- escapeBraces(paste0(skipped, " (", reasons, ")", collapse = ", "))
            cli::cli_warn(c(
                "!" = "{.field {i}}: no threshold for {length(skipped)} of
                       {length(levs)} {.field {group_by}} {cli::qty(length(levs))}group{?s}.",
                "i" = "Skipped: {detail}",
                "x" = "Cells in {cli::qty(length(skipped))}{?that group/those groups}
                       are left as {.val {NA}}, not {.val {FALSE}}."
            ))
        }

        reportFlagMetric(i, metric_df, use_lo, use_hi, is.null(group_by), n_cells)

        if (any(!is.na(metric_df$column))) {
            obj[[column]] <- stats::setNames(outliers, colnames(obj))
        }
        rows[[length(rows) + 1L]] <- metric_df
    }

    summary_df <- do.call(rbind, rows)
    rownames(summary_df) <- NULL

    #== Keep the table describing every flag column on the object, not just the
    #== ones this call happened to write.
    prev <- obj@misc[["flagQC"]]
    if (is.data.frame(prev) && identical(names(prev), names(summary_df))) {
        prev <- prev[!prev$column %in% stats::na.omit(summary_df$column), , drop = FALSE]
        summary_df <- rbind(prev, summary_df)
        rownames(summary_df) <- NULL
    }
    obj@misc[["flagQC"]] <- summary_df

    written <- length(unique(stats::na.omit(summary_df$column)))
    cli::cli_alert_success(
        "Added {written} outlier column{?s}; summary in {.code {arg_name}@misc$flagQC}."
    )

    obj
}

#' Retrieve the Outlier Summary Table from an Object
#'
#' Returns the table [flagQC()] recorded when it flagged the object: one row
#' per metric per group, holding the thresholds that were applied and how many
#' cells each one caught. It is a plain data frame, so it can be sorted,
#' subset, written to CSV or joined against a sample sheet like any other.
#'
#' The table is a record of what was flagged, not a live view of the object. It
#' is written when [flagQC()] runs and is not recalculated afterwards, so cell
#' counts describe the object as it stood at that moment. Subsetting the object
#' does not update them -- which is the point, since the whole reason to keep
#' the table is to be able to say afterwards how many cells were dropped and
#' why.
#'
#' Repeated [flagQC()] calls accumulate into one table, so an object flagged on
#' counts globally and then on mitochondrial content per sample yields a single
#' table describing both. Rows for a flag column that was overwritten are
#' replaced rather than duplicated, so every row corresponds to a metadata
#' column that is actually present.
#'
#' @param obj A [Seurat][Seurat::CreateSeuratObject] object that has been
#'   through [flagQC()].
#'
#' @return A data frame with one row per metric per group. See the *Value*
#'   section of [flagQC()] for what each column means. Aborts if the object has
#'   never been flagged, rather than returning `NULL` for a caller to trip
#'   over.
#'
#' @examples
#' \dontrun{
#' seurat <- diagnoseQC(seurat, organism = "mmu")
#' seurat <- flagQC(
#'     seurat,
#'     metrics  = c(nCount_RNA = "lower", pctMito = "higher"),
#'     group_by = "Ident"
#' )
#'
#' qc <- getOutlierSummary(seurat)
#' qc
#'
#' # the thresholds that were actually applied, per sample
#' qc[, c("metric", "group", "lower", "upper", "n_lower", "n_upper")]
#'
#' # which samples lose the most cells, and to which metric?
#' qc[order(-qc$pct_flagged), c("metric", "group", "n_flagged", "pct_flagged")]
#'
#' # groups no threshold could be estimated for
#' qc[is.na(qc$column), c("metric", "group", "n_cells")]
#'
#' # keep it alongside the filtered object
#' write.csv(qc, "qc/outlier_summary.csv", row.names = FALSE)
#' }
#'
#' @seealso [flagQC()], which writes this table, and [diagnoseQC()], which
#'   calculates the metrics it describes.
#'
#' @export
getOutlierSummary <- function(obj) {
    arg_name <- deparse(substitute(obj))
    validateQCObject(obj, arg_name, fn = "getOutlierSummary")

    summary_df <- obj@misc[["flagQC"]]

    if (is.null(summary_df)) {
        cli::cli_abort(c(
            "{.arg {arg_name}} has no outlier summary.",
            "i" = "{.fn flagQC} writes one to {.code @misc$flagQC}; run it first."
        ))
    }

    #== Anything else under that name was not written by flagQC(), and handing
    #== it back as though it were the summary table would be worse than saying so.
    if (!is.data.frame(summary_df)) {
        cli::cli_abort(c(
            "{.code {arg_name}@misc$flagQC} is {.cls {class(summary_df)[1]}},
             not a data frame.",
            "i" = "It was overwritten by something other than {.fn flagQC}."
        ))
    }

    summary_df
}

#' Resolve `group_by` to a factor of one level per cell
#'
#' Returns a single-level factor when grouping is off, so the caller loops over
#' groups either way instead of branching on a mode flag.
#'
#' @param group_by The caller's argument, unmodified.
#' @param md The object's `meta.data`.
#' @param arg_name Name the object was passed under, for the message.
#' @param n_cells Cells in the object.
#' @param call Environment to blame in the error.
#'
#' @return A factor of length `n_cells`, possibly containing `NA`.
#'
#' @noRd
resolveFlagGroups <- function(group_by, md, arg_name, n_cells,
                              call = parent.frame()) {
    if (is.null(group_by)) {
        return(factor(rep("all", n_cells)))
    }
    if (!is.character(group_by) || length(group_by) != 1 || is.na(group_by)) {
        cli::cli_abort(
            "{.arg group_by} must be a single metadata column name, or {.code NULL}.",
            call = call
        )
    }
    if (!group_by %in% colnames(md)) {
        cli::cli_abort(c(
            "{.arg group_by} column {.val {group_by}} is not in {.arg {arg_name}}.",
            "i" = "{.field Ident} is set to the sample name by {.fn seuratFrom10x}."
        ), call = call)
    }

    g <- md[[group_by]]

    #== A continuous column would give one group per cell and no thresholds at
    #== all. Integer codes are a legitimate way to label samples, so only
    #== doubles are refused.
    if (is.double(g)) {
        cli::cli_abort(c(
            "{.arg group_by} column {.val {group_by}} is continuous
             ({.cls {class(g)[1]}}), so it has no groups to estimate within.",
            "i" = "Group on a discrete column, or convert this one with
                   {.fn factor} first."
        ), call = call)
    }

    groups <- factor(as.character(g))
    if (length(levels(groups)) == 0) {
        cli::cli_abort(
            "{.arg group_by} column {.val {group_by}} is entirely missing.",
            call = call
        )
    }

    n_na <- sum(is.na(groups))
    if (n_na > 0) {
        cli::cli_warn(c(
            "!" = "{n_na} cell{?s} {?has/have} no {.field {group_by}} value.",
            "i" = "{cli::qty(n_na)}{?It/They} belong{?s/} to no group, so
                   {cli::qty(n_na)}{?its/their} flag{?s} {?is/are} left as
                   {.val {NA}}."
        ))
    }

    groups
}

#' Median and MAD thresholds for one vector, on the vector's own scale
#'
#' Estimation happens on the transformed scale; the thresholds are transformed
#' back before being returned, so the caller compares raw values against them
#' and the numbers it reports are the numbers it applied.
#'
#' @param x Numeric vector, possibly containing `NA`.
#' @param nmads Number of MADs from the median.
#' @param log Estimate on a `log2` scale.
#'
#' @return A list with `ok = FALSE` and a `reason`, or `ok = TRUE` plus
#'   `center`, `lower` and `upper` on the input's own scale.
#'
#' @noRd
madThresholds <- function(x, nmads, log) {
    vals <- x[!is.na(x)]

    #== mad() of two values is arithmetic on a coin flip, and of one value is
    #== always zero. Three is the floor at which the statistic means anything.
    if (length(vals) < 3) {
        return(list(ok = FALSE, reason = "fewer than 3 cells"))
    }

    shift <- if (log && any(vals == 0)) 1 else 0
    xt    <- if (log) log2(vals + shift) else vals

    med <- stats::median(xt)
    dev <- stats::mad(xt, center = med)

    #== A zero MAD means over half the cells share one value, so both
    #== thresholds collapse onto it and every flag would be an artefact.
    if (!is.finite(dev) || dev == 0) {
        return(list(ok = FALSE, reason = "zero MAD"))
    }

    back <- if (log) function(v) 2^v - shift else identity
    list(ok     = TRUE,
         reason = NA_character_,
         center = back(med),
         lower  = back(med - nmads * dev),
         upper  = back(med + nmads * dev))
}

#' Report one metric's flagging to the console
#'
#' Grouped runs roll the per-group rows up into a single line plus one bullet
#' per bound, showing the range of thresholds across groups. Printing every
#' group would put a 20-sample object's console output past a screenful per
#' metric, and the per-group detail is in the returned table anyway.
#'
#' @param metric Metric name.
#' @param df The metric's rows of the summary table.
#' @param use_lo,use_hi Whether each bound was applied.
#' @param global Whether thresholds were estimated globally.
#' @param n_cells Cells in the object.
#'
#' @return `invisible(NULL)`, called for the console output.
#'
#' @noRd
reportFlagMetric <- function(metric, df, use_lo, use_hi, global, n_cells) {
    n_flag <- sum(df$n_flagged, na.rm = TRUE)
    n_used <- sum(df$n_cells[!is.na(df$column)])
    pct    <- if (n_used > 0) round(100 * n_flag / n_used, 2) else NA_real_

    if (n_used == 0) {
        cli::cli_alert_info(
            "{.field {metric}}: no cells judged -- see the warning above."
        )
        return(invisible(NULL))
    }

    scope <- if (global) {
        ""
    } else {
        cli::format_inline(" across {sum(!is.na(df$column))} group{?s}")
    }
    cli::cli_alert_info(
        "{.field {metric}}: {n_flag}/{n_used} cell{?s} flagged ({pct}%){scope}."
    )

    #== One threshold reads as a number; many read as the range they span,
    #== which is the thing worth seeing -- groups needing wildly different
    #== cutoffs are the reason to be grouping in the first place.
    fmt <- function(v) {
        v <- v[!is.na(v)]
        if (length(v) == 0) return("NA")
        rng <- signif(range(v), 4)
        if (rng[1] == rng[2]) as.character(rng[1]) else paste0(rng[1], "-", rng[2])
    }

    bullets <- character(0)
    if (use_lo) {
        bullets <- c(bullets, "*" = cli::format_inline(
            "lower: < {fmt(df$lower)} -- {sum(df$n_lower, na.rm = TRUE)} cell{?s}"
        ))
    }
    if (use_hi) {
        bullets <- c(bullets, "*" = cli::format_inline(
            "upper: > {fmt(df$upper)} -- {sum(df$n_upper, na.rm = TRUE)} cell{?s}"
        ))
    }
    n_missing <- sum(df$n_missing, na.rm = TRUE)
    if (n_missing > 0) {
        bullets <- c(bullets, "!" = cli::format_inline(
            "{n_missing} cell{?s} with a missing value, left as {.val {NA}}"
        ))
    }
    cli::cli_bullets(escapeBraces(bullets))

    invisible(NULL)
}

#' One row of the flagQC() summary table
#'
#' Built in a helper so the skipped-group row and the ordinary row cannot drift
#' apart in column name, order or type -- `rbind()` on data frames whose
#' columns disagree is either an error or, worse, a silent coercion.
#'
#' @noRd
flagQCRow <- function(metric, bound, nmads, log, group_by, group,
                      center, lower, upper, n_lower, n_upper, n_flagged,
                      n_missing, n_cells, column) {
    data.frame(
        metric      = metric,
        bound       = bound,
        nmads       = as.numeric(nmads),
        scale       = if (log) "log2" else "identity",
        group_by    = if (is.null(group_by)) NA_character_ else group_by,
        group       = if (is.null(group_by)) NA_character_ else as.character(group),
        center      = as.numeric(center),
        lower       = as.numeric(lower),
        upper       = as.numeric(upper),
        n_lower     = as.integer(n_lower),
        n_upper     = as.integer(n_upper),
        n_flagged   = as.integer(n_flagged),
        pct_flagged = if (n_cells > 0) 100 * as.integer(n_flagged) / n_cells else NA_real_,
        n_missing   = as.integer(n_missing),
        n_cells     = as.integer(n_cells),
        column      = column,
        stringsAsFactors = FALSE
    )
}
