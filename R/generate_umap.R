#' Generate UMAP plot
#' 
#' Generates the matrix of UMAP visualisation coordinates from a given 
#' SummarizedExperiment.
#' The modified SummarizedExperiment object is returned, with UMAP coordinates
#' stored in the `reducedDims()` slot.
#' 
#' @param summarizedExperiment the SingleCellExperiment of SummarizedExperiment
#' object containing the counts to plot
#' @param counts character string indicating the name of the assay to generate UMAP visualtization from.
#' The assay should be accessable by `assay(summarizedExperiemnt, counts)`.
#' @seealso [sc_long_pipeline()]
#' @seealso [bulk_long_pipeline()]
#' @importFrom scater runUMAP
#' @importFrom SummarizedExperiment assayNames
#' @importFrom dplyr expr
#' 
#' @return Modified SummarizedExperiment (or SingleCellExperiment) containing
#' the UMAP visualisation coordinates in the `reducedDims()` UMAP slot.
generate_umap <- function(summarizedExperiment, counts="counts") {
  if (!(counts %in% assayNames(summarizedExperiment))) {
    stop(paste0(counts, " not found in assays"))
  } 
  summarizedExperiment <- runUMAP(summarizedExperiment, expr)
  
  return(summarizedExperiment)
}