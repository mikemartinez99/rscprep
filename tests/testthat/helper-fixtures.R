# Test fixtures for rscprep.
#
# Strategy: the package's real logic is *path resolution* -- deciding which H5
# file to look for, what to name the resulting sample, and whether a missing
# file is fatal. None of that requires a real H5 file, so fixtures create the
# Cell Ranger directory layouts with empty placeholder files (the code only
# ever calls file.exists() on them) and mock Seurat::Read10X_h5() to return an
# in-memory counts matrix. Seurat::CreateSeuratObject() is NOT mocked -- it
# runs for real, so the tests also verify we hand it valid arguments.

# ---- counts matrix -------------------------------------------------------

#' Build a small sparse counts matrix.
#'
#' With `extras = TRUE` the matrix contains two deliberately marginal entries
#' used to test min.cells / min.features pass-through:
#'   * `RareGene`   -- detected in only 2 cells (dropped by min.cells >= 3)
#'   * `SparseCell` -- only 5 features detected (dropped by min.features >= 6)
#' Everything else is dense enough to survive the package defaults
#' (min.cells = 10, min.features = 50), so the default-argument case yields a
#' predictable 60 features x 20 cells.
make_counts <- function(n_genes = 60, n_cells = 20, extras = TRUE) {
    m <- matrix(
        1,
        nrow     = n_genes,
        ncol     = n_cells,
        dimnames = list(paste0("Gene", seq_len(n_genes)),
                        paste0("Cell", seq_len(n_cells)))
    )

    if (extras) {
        rare <- matrix(0, nrow = 1, ncol = ncol(m),
                       dimnames = list("RareGene", colnames(m)))
        rare[1, 1:2] <- 1
        m <- rbind(m, rare)

        sparse <- matrix(0, nrow = nrow(m), ncol = 1,
                         dimnames = list(rownames(m), "SparseCell"))
        sparse[1:5, 1] <- 1
        m <- cbind(m, sparse)
    }

    Matrix::Matrix(m, sparse = TRUE)
}

# ---- directory layouts ---------------------------------------------------

#' Create a `cellranger count` layout: <root>/<sample>/outs/filtered_...h5
#' Returns the `outs` directory path.
make_count_run <- function(root, sample = "Sample_1", h5 = TRUE) {
    outs <- file.path(root, sample, "outs")
    dir.create(outs, recursive = TRUE, showWarnings = FALSE)
    if (h5) file.create(file.path(outs, "filtered_feature_bc_matrix.h5"))
    outs
}

#' Create a `cellranger multi` layout:
#'   <root>/<run>/outs/per_sample_outs/<sample>/count/sample_filtered_...h5
#'
#' `where = "count"` is the real Cell Ranger layout. `where = "flat"` drops the
#' H5 straight into the sample directory instead, which is what staging or
#' flattening pipelines produce and what the fallback candidate exists for.
#'
#' `h5` may be a single logical or one logical per sample, so a run can be
#' built with some samples' H5 files deliberately absent.
#' Returns the `outs` directory path.
make_multi_run <- function(root, run = "Pool1", samples = c("Sample_A", "Sample_B"),
                           h5 = TRUE, make_per_sample_outs = TRUE,
                           where = c("count", "flat")) {
    where <- match.arg(where)
    outs  <- file.path(root, run, "outs")
    dir.create(outs, recursive = TRUE, showWarnings = FALSE)

    if (make_per_sample_outs) {
        h5 <- rep_len(h5, length(samples))
        for (i in seq_along(samples)) {
            sdir <- file.path(outs, "per_sample_outs", samples[i])
            dir.create(sdir, recursive = TRUE, showWarnings = FALSE)
            if (h5[i]) {
                file.create(sampleFilteredH5(sdir, where))
            }
        }
    }

    outs
}

#' Path of a demuxed sample's filtered H5 within its per_sample_outs directory.
#'
#' Creates the `count/` folder as a side effect when `where = "count"`, so
#' callers can hand the result straight to file.create().
sampleFilteredH5 <- function(sample_dir, where = c("count", "flat")) {
    where <- match.arg(where)

    dir <- if (where == "count") file.path(sample_dir, "count") else sample_dir
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)

    file.path(dir, "sample_filtered_feature_bc_matrix.h5")
}

#' Add a raw feature barcode matrix to an existing `outs` directory.
#'
#' `where = "flat"` writes `<outs>/raw_feature_bc_matrix.h5`, the
#' `cellranger count` location. `where = "multi"` writes
#' `<outs>/multi/count/raw_feature_bc_matrix.h5`, where `cellranger multi`
#' puts the single run-level raw matrix it produces for the whole pool.
#' Returns the path of the file created.
add_raw_h5 <- function(outs, where = c("flat", "multi")) {
    where <- match.arg(where)

    dir <- if (where == "multi") file.path(outs, "multi", "count") else outs
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)

    h5 <- file.path(dir, "raw_feature_bc_matrix.h5")
    file.create(h5)
    h5
}

#' A raw-matrix-shaped counts matrix with a knee in it.
#'
#' `n_cells` barcodes at ~`depth` UMIs sit above a long tail of `n_empty`
#' ambient barcodes, which is the shape DropletUtils::barcodeRanks() needs to
#' find a knee and an inflection point. Kept small enough to rank quickly, and
#' deliberately spread over many distinct totals so the fit has unique points
#' to work with.
make_raw_counts <- function(n_genes = 40, n_cells = 300, n_empty = 3000,
                            depth = 40, seed = 42) {
    withr::with_seed(seed, {
        cells <- Matrix::rsparsematrix(
            n_genes, n_cells, density = 0.7,
            rand.x = function(n) stats::rpois(n, depth)
        )
        empty <- Matrix::rsparsematrix(
            n_genes, n_empty, density = 0.05,
            rand.x = function(n) stats::rpois(n, 1) + 1
        )
    })

    m <- cbind(cells, empty)
    dimnames(m) <- list(
        paste0("Gene", seq_len(n_genes)),
        paste0(sprintf("%016d", seq_len(n_cells + n_empty)), "-1")
    )
    methods::as(m, "CsparseMatrix")
}

# ---- mocking -------------------------------------------------------------

#' Mock Seurat::Read10X_h5() for the duration of the calling test.
#'
#' Returns an environment whose `$files` field accumulates every `filename`
#' the package asked for, in call order -- this is how tests assert that the
#' correct H5 path was resolved, not merely that *a* file was read.
local_mock_read10x <- function(counts = make_counts(), .env = parent.frame()) {
    calls <- new.env(parent = emptyenv())
    calls$files <- character(0)

    testthat::local_mocked_bindings(
        Read10X_h5 = function(filename, ...) {
            calls$files <- c(calls$files, filename)
            counts
        },
        .package = "Seurat",
        .env     = .env
    )

    calls
}

#' A temp directory scoped to the calling test.
local_root <- function(.env = parent.frame()) {
    withr::local_tempdir(.local_envir = .env)
}

# ---- in-memory Seurat objects --------------------------------------------
#
# mergeSeurat() takes objects, not paths, so its fixtures skip the filesystem
# entirely and build real Seurat objects directly.

#' Build a Seurat object with a known shape.
#'
#' `ident` sets both `project` and the `Ident` metadata column, mirroring what
#' seuratFrom10x() produces. `ident = NULL` omits `Ident` altogether, which is
#' how the post-merge Ident check is exercised.
make_seurat <- function(ident = "Sample_A", n_cells = 10, n_genes = 20,
                        assay = "RNA", cells = NULL) {
    if (is.null(cells)) cells <- paste0("Cell", seq_len(n_cells))

    m <- matrix(
        1,
        nrow     = n_genes,
        ncol     = length(cells),
        dimnames = list(paste0("Gene", seq_len(n_genes)), cells)
    )

    obj <- Seurat::CreateSeuratObject(
        counts       = Matrix::Matrix(m, sparse = TRUE),
        project      = if (is.null(ident)) "SeuratProject" else ident,
        min.cells    = 0,
        min.features = 0,
        assay        = assay
    )

    if (!is.null(ident)) obj$Ident <- ident
    obj
}

#' A named list of Seurat objects, one per ident.
make_seurat_list <- function(idents = c("Sample_A", "Sample_B"), ...) {
    stats::setNames(lapply(idents, make_seurat, ...), idents)
}

#' A Seurat object with zero cells.
#'
#' Seurat refuses to *construct* one (CreateSeuratObject and `[` both reject
#' empty input), so this hollows out a real object instead: emptied count and
#' data matrices, no metadata rows. The result still passes `inherits()` and
#' `validObject()`, which is exactly what mergeSeurat()'s zero-cell guard is
#' there to catch.
make_empty_seurat <- function(ident = "Sample_Empty", n_genes = 20) {
    obj   <- make_seurat(ident, n_cells = 2, n_genes = n_genes)
    empty <- Matrix::Matrix(
        matrix(0, nrow = n_genes, ncol = 0,
               dimnames = list(rownames(obj), character(0))),
        sparse = TRUE
    )

    obj@assays[[Seurat::DefaultAssay(obj)]]@counts <- empty
    obj@assays[[Seurat::DefaultAssay(obj)]]@data   <- empty
    obj@meta.data    <- obj@meta.data[0, , drop = FALSE]
    obj@active.ident <- factor(character(0))

    obj
}

#' Build a sample-level metadata sheet.
#'
#' `key = "column"` puts sample names in a `Sample` column (how a CSV round
#' trips); `key = "rownames"` puts them in the rownames. `...` supplies the
#' metadata columns, recycled to one value per sample.
make_sample_sheet <- function(idents = c("Sample_A", "Sample_B"),
                              key = c("column", "rownames"), ...) {
    key  <- match.arg(key)
    cols <- list(...)
    if (length(cols) == 0) cols <- list(Treatment = c("Control", "Bleomycin"))

    cols <- lapply(cols, rep_len, length(idents))
    df   <- as.data.frame(cols, stringsAsFactors = FALSE)

    if (key == "column") {
        df <- cbind(Sample = idents, df, stringsAsFactors = FALSE)
    } else {
        rownames(df) <- idents
    }

    df
}

#' Something that claims to be a Seurat object but is not a valid one.
#'
#' `inherits(x, "Seurat")` is TRUE, `validObject(x)` fails -- the case the
#' two-part validity check in mergeSeurat() exists for.
make_fake_seurat <- function() {
    structure(list(), class = "Seurat")
}

# ---- QC fixtures ---------------------------------------------------------
#
# diagnoseQC() and flagQC() work on gene *names* and metadata *values*, not on
# file layouts, so their fixtures build objects with both pinned exactly.
#
# Every count in make_pattern_seurat() is 1, which makes the percentage
# metrics exact rather than approximate: with one count per gene per cell,
# pctMito is 100 * n_mito / n_genes for every cell, so the expected value is
# arithmetic rather than a tolerance.

#' A Seurat object whose gene names carry organism-specific QC patterns.
#'
#' Defaults are mouse-style symbols (`mt-`, `Rp[sl]`). Pass human-style names
#' to exercise a different organism, or `character(0)` to build an object a
#' pattern cannot match at all. `extra` adds further named genes -- cell cycle
#' markers, say -- alongside the anonymous background.
make_pattern_seurat <- function(mito    = c("mt-Co1", "mt-Co2", "mt-Nd1"),
                                ribo    = c("Rps3", "Rpl4"),
                                n_other = 15,
                                n_cells = 20,
                                ident   = "Sample_A",
                                assay   = "RNA",
                                extra   = character(0)) {
    genes <- c(mito, ribo, extra, paste0("Gene", seq_len(n_other)))
    cells <- paste0("Cell", seq_len(n_cells))

    m <- matrix(1, nrow = length(genes), ncol = n_cells,
                dimnames = list(genes, cells))

    obj <- Seurat::CreateSeuratObject(
        counts       = Matrix::Matrix(m, sparse = TRUE),
        project      = ident,
        min.cells    = 0,
        min.features = 0,
        assay        = assay
    )
    obj$Ident <- ident
    obj
}

#' A Seurat object carrying metadata columns with exactly these values.
#'
#' flagQC() reads metadata, never counts, so setting the columns directly is
#' what makes its thresholds predictable: the median and MAD of the input are
#' known, so the expected cutoffs are arithmetic.
#'
#' `...` names the metadata columns; every one must be the same length, which
#' sets the object's cell count. `ident` is recycled to one value per cell, so
#' a vector splits the object into groups for `group_by`.
make_metric_seurat <- function(..., ident = "Sample_A") {
    metrics <- list(...)
    if (length(metrics) == 0) stop("supply at least one metric column")

    n   <- length(metrics[[1]])
    obj <- make_seurat(ident = "Sample_A", n_cells = n, n_genes = 10)

    for (nm in names(metrics)) {
        obj[[nm]] <- stats::setNames(metrics[[nm]], colnames(obj))
    }
    obj$Ident <- stats::setNames(rep_len(ident, n), colnames(obj))
    obj
}

#' Twenty values with a median of 6 and a MAD of exactly 2 * 1.4826.
#'
#' Sorted absolute deviations from the median are 0,0,0,0 then 2 x8 then 4 x8,
#' so the median deviation is 2 whether or not `extreme` is appended -- adding
#' one large value shifts neither statistic, which is the property that makes
#' the expected thresholds hand-computable.
#'
#' With `nmads = 3` and `log = FALSE` the thresholds are 6 +/- 8.8956, so
#' `extreme = 100` is flagged high and nothing is flagged low.
metric_values <- function(extreme = NULL) {
    c(rep(c(2, 4, 6, 8, 10), each = 4), extreme)
}

#' The same shape on a log2 scale: values are powers of two, so log2() of them
#' is 1,2,3,4,5 (four each) and the median and MAD are exact there too.
#'
#' log2 median is 3 and log2 MAD is 1.4826, so with `nmads = 3` the upper
#' threshold is 2^(3 + 3 * 1.4826) ~ 174.7 and `extreme = 1024` is flagged.
log_metric_values <- function(extreme = NULL) {
    c(rep(c(2, 4, 8, 16, 32), each = 4), extreme)
}

#' The MAD constant stats::mad() scales by.
MAD_CONST <- 1.4826

# ---- cell cycle fixtures -------------------------------------------------
#
# Seurat::CellCycleScoring() bins genes by expression and samples 100 control
# genes per bin, so it needs a background pool of at least ~1,200 genes before
# it will run at all -- below that it fails with "cannot take a sample larger
# than the population", which is the error runCellCycle() wraps. The fixture
# is sized just above that floor and kept narrow in cells to stay quick.

#' A Seurat object the cell cycle markers can actually be scored against.
#'
#' `markers` supplies the marker symbols to embed; the default is the human
#' set Seurat ships, so `organism = "hsa"` scores without any ortholog
#' conversion. Pass converted symbols to build a non-human object, or a subset
#' to exercise the few-markers-matched path.
#'
#' `normalise = FALSE` returns raw counts, which is what the
#' normalise-on-the-way-through behaviour is tested against.
make_cellcycle_seurat <- function(markers   = c(Seurat::cc.genes$s.genes,
                                                Seurat::cc.genes$g2m.genes),
                                  n_cells   = 40,
                                  n_bg      = 1200,
                                  normalise = TRUE,
                                  seed      = 42) {
    genes <- c(markers, paste0("BG", seq_len(n_bg)))

    m <- withr::with_seed(seed, {
        matrix(stats::rpois(length(genes) * n_cells, 5),
               nrow     = length(genes),
               dimnames = list(genes, paste0("Cell", seq_len(n_cells))))
    })

    obj <- Seurat::CreateSeuratObject(Matrix::Matrix(m, sparse = TRUE),
                                      min.cells = 0, min.features = 0)
    if (normalise) obj <- Seurat::NormalizeData(obj, verbose = FALSE)
    obj
}

#' Mock orthogene::convert_orthologs() for the duration of the calling test.
#'
#' Maps every input symbol to a predictable pseudo-ortholog (`MCM5` becomes
#' `Xmcm5`) so the conversion is deterministic and offline. `drop` names input
#' symbols to omit from the result, standing in for genes with no one-to-one
#' ortholog -- which is what makes the output a different length and order
#' from the input, and therefore what the match-by-symbol logic is for.
#'
#' Returns an environment whose `$calls` field records each `gene_df` handed
#' over, so a test can assert both sets were converted in a single request.
local_mock_orthologs <- function(drop = character(0), shuffle = TRUE,
                                 .env = parent.frame()) {
    calls <- new.env(parent = emptyenv())
    calls$gene_df <- list()
    calls$output_species <- character(0)

    testthat::local_mocked_bindings(
        convert_orthologs = function(gene_df, output_species, ...) {
            calls$gene_df <- c(calls$gene_df, list(gene_df))
            calls$output_species <- c(calls$output_species, output_species)

            kept <- setdiff(gene_df, drop)
            if (shuffle) kept <- rev(kept)   # order must not be relied on

            data.frame(input_gene    = kept,
                       ortholog_gene = paste0("X", tolower(kept)),
                       stringsAsFactors = FALSE)
        },
        .package = "orthogene",
        .env     = .env
    )

    calls
}

#' The pseudo-ortholog symbol local_mock_orthologs() maps a human symbol to.
pseudo_ortholog <- function(x) paste0("X", tolower(x))
