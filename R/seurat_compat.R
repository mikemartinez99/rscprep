# Seurat 4 / Seurat 5 compatibility
#
# Every version difference the package has to handle goes through here, so a
# function never branches on the Seurat version itself. Three differences
# matter:
#
#   * Addressing assay data. SeuratObject 4 takes `slot = "data"`; 5 renamed it
#     `layer`, and has since made `slot` defunct, so the old name errors.
#   * What an unnormalised object holds. SeuratObject 4 copies counts into the
#     `data` slot at creation. A 5 (Assay5) object has no `data` layer at all
#     until NormalizeData() runs.
#   * What merge() returns. SeuratObject 4 binds counts into one matrix. 5
#     keeps one layer per input (`counts.Sample_A`, `counts.Sample_B`, ...).
#
# The version is read from SeuratObject rather than Seurat because the classes
# and accessors live there. An object's own class can still differ from the
# installed version -- an object saved under Seurat 4 loads under 5 as a
# classic `Assay`, not an `Assay5` -- so the helpers also check the assay
# itself where that matters.
#
# Functions that exist only in SeuratObject 5 are fetched with
# getExportedValue() rather than `::`, so the package still installs and
# passes R CMD check against SeuratObject 4.

#' Is SeuratObject 5 or later installed?
#'
#' @return `TRUE` or `FALSE`.
#'
#' @noRd
isSeuratV5 <- function() {
    utils::packageVersion("SeuratObject") >= "5.0.0"
}

#' A SeuratObject 5 function, fetched at call time
#'
#' @noRd
seuratV5Fn <- function(name) {
    getExportedValue("SeuratObject", name)
}

#' Read one kind of assay data, whichever Seurat is installed
#'
#' @param obj A Seurat object.
#' @param layer `"counts"`, `"data"` or `"scale.data"`.
#' @param assay Assay to read. Defaults to the object's default assay.
#'
#' @return A matrix with features as rows and every cell as columns, or `NULL`
#'   when the assay holds no such data (only possible under Seurat 5, e.g. the
#'   `data` layer of an object that has not been normalised). Layers split by
#'   sample after a Seurat 5 merge are bound back together, so the result
#'   always covers every cell.
#'
#' @noRd
getLayerValues <- function(obj, layer, assay = NULL) {
    if (is.null(assay)) assay <- Seurat::DefaultAssay(obj)

    if (!isSeuratV5()) {
        return(Seurat::GetAssayData(obj, slot = layer, assay = assay))
    }

    layers <- seuratV5Fn("Layers")(obj[[assay]])
    if (layer %in% layers) {
        return(seuratV5Fn("LayerData")(obj[[assay]], layer = layer))
    }

    split <- layers[startsWith(layers, paste0(layer, "."))]
    if (length(split) == 0) return(NULL)

    values <- lapply(split, function(l) seuratV5Fn("LayerData")(obj[[assay]], layer = l))
    do.call(cbind, values)
}

#' Join split layers back into one per kind of data
#'
#' A no-op under Seurat 4, whose merge() never splits, and for any assay that
#' is not an `Assay5` (a classic assay has no layers to join).
#'
#' @param obj A Seurat object.
#'
#' @return `obj`, with every `Assay5` assay's layers joined.
#'
#' @noRd
joinAssayLayers <- function(obj) {
    if (!isSeuratV5()) return(obj)

    for (assay in names(obj@assays)) {
        if (methods::is(obj[[assay]], "Assay5")) {
            obj[[assay]] <- seuratV5Fn("JoinLayers")(obj[[assay]])
        }
    }
    obj
}
