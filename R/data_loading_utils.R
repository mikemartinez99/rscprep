#' @noRd
readOneSample10x <- function(outs_dir, ident, min.cells, min.features,
                              filename = "filtered_feature_bc_matrix.h5",
                              on_missing = warning) {
    mtx_path <- file.path(outs_dir, filename)
    if (!file.exists(mtx_path)) {
        on_missing(sprintf("File not found at '%s'", mtx_path))
        return(NULL)
    }
    counts <- Seurat::Read10X_h5(filename = mtx_path)
    seurat <- Seurat::CreateSeuratObject(
        counts       = counts,
        project      = ident,
        min.cells    = min.cells,
        min.features = min.features,
        assay        = "RNA"
    )
    seurat$Ident <- ident
    seurat
}

#' @noRd
readMultiplexedRun10x <- function(outs_dir, min.cells, min.features, on_missing = warning) {
    per_sample_dir <- file.path(outs_dir, "per_sample_outs")
    individual_samples <- list.dirs(per_sample_dir, full.names = TRUE, recursive = FALSE)

    sample_seurats <- lapply(individual_samples, function(i) {
        ident <- basename(i)
        readOneSample10x(
            outs_dir     = i,
            ident        = ident,
            min.cells    = min.cells,
            min.features = min.features,
            filename     = "sample_filtered_feature_bc_matrix.h5",
            on_missing   = on_missing
        )
    })
    names(sample_seurats) <- basename(individual_samples)
    sample_seurats
}
