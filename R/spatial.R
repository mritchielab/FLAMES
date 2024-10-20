#' Create a SpatialExperiment object
#'
#' @description This function creates a SpatialExperiment object from a SingleCellExperiment object and a spatial barcode file.
#' @param sce The SingleCellExperiment object obtained from running the \code{\link{sc_long_pipeline}} function.
#' @param spatial_barcode_file The path to the spatial barcode file, e.g. \code{"spaceranger-2.1.1/lib/python/cellranger/barcodes/visium-v2_coordinates.txt"}.
#' @param mannual_align_json The path to the mannual alignment json file.
#' @param image 'DataFrame' containing the image data. See \code{?SpatialExperiment::readImgData} and \code{?SpatialExperiment::SpatialExperiment}.
#' @return A SpatialExperiment object.
#' @importFrom readr read_table
#' @importFrom jsonlite fromJSON
#' @importFrom tidyr as_tibble
#' @importFrom dplyr mutate left_join select
#' @importFrom SummarizedExperiment assays colData rowData rowRanges rowRanges<- colData<-
#' @importFrom SpatialExperiment SpatialExperiment readImgData imgData imgData<-
#' @export
create_spe <- function(sce, spatial_barcode_file, mannual_align_json, image) {

  # Read the full list file
  full_list <- readr::read_table(
    spatial_barcode_file,
    col_names = c("barcode", "col", "row"),
    col_types = "cii"
  )

  # use mannual alignment of image and spots if provided
  if (!missing(mannual_align_json)) {
    align_df <- jsonlite::fromJSON(mannual_align_json)$oligo |>
      tidyr::as_tibble() |>
      dplyr::mutate(row = row + 1, col = col + 1)
    full_list <- dplyr::left_join(align_df, full_list, by = c('row', 'col'))
  }

  # add spatial info to colData
  coldata <- full_list |>
    dplyr::select(-barcode) |>
    as.data.frame()
  rownames(coldata) <- full_list$barcode
  coldata <- coldata[colnames(sce), ]
  coldata <- cbind(SummarizedExperiment::colData(sce), coldata)

  # Create a SpatialExperiment object
  spe <- SpatialExperiment::SpatialExperiment(
    assay = SummarizedExperiment::assays(sce),
    colData = coldata,
    rowData = SummarizedExperiment::rowData(sce),
    spatialCoordsNames = c("imageX", "imageY")
  )
  SummarizedExperiment::rowRanges(spe) <- SummarizedExperiment::rowRanges(sce)

  if ("tissue" %in% names(SummarizedExperiment::colData(spe))) {
    SummarizedExperiment::colData(spe)$in_tissue <- SummarizedExperiment::colData(spe)$tissue
    cat("in_tissue column is set to tissue column from", mannual_align_json, "\n")
    cat("in_tissue counts:")
    print(table(SummarizedExperiment::colData(spe)$in_tissue, useNA = "always"))
  } else {
    SummarizedExperiment::colData(spe)$in_tissue <- TRUE
    warning("No tissue column in full list file, all spots are considered to be in tissue")
  }

  # add image file if provided
  if (!missing(image)) {
    if (is.character(image)) {
      image <- SpatialExperiment::readImgData(
        image,
        sample_id = SummarizedExperiment::colData(spe)$sample_id[1],
        imageSources = file.path(image, "tissue_hires_image.png")
      )
    }
    SpatialExperiment::imgData(spe) <- image
  }

  return(spe)
}

#' Plot spatial pie chart
#'
#' @description This function plots a spatial pie chart for given features.
#' @param spe The SpatialExperiment object.
#' @param features The features to plot.
#' @param assay_type The assay that contains the given features.
#' @importFrom ggplot2 ggplot annotation_raster coord_fixed theme_void aes
#' @importFrom scatterpie geom_scatterpie
#' @importFrom dplyr mutate
#' @return A ggplot object.
plot_spatial_pie <- function(spe, features, assay_type = 'counts') {
  if (nrow(imgData(spe)) > 0) {
    background_img <- SpatialExperiment::imgData(spe)$data[[1]] |>
      SpatialExperiment::imgRaster()
    maxX <- dim(background_img)[1]
    maxY <- dim(background_img)[2]
    p1 <- ggplot2::ggplot(mapping = ggplot2::aes(1:maxX, 1:maxY)) +
      ggplot2::annotation_raster(background_img, 
        xmin = 1, xmax = maxX, ymin = 1, ymax = maxY)
  } else {
    maxX <- max(plot_d$imageX)
    minX <- min(plot_d$imageX)
    maxY <- max(plot_d$imageY)
    minY <- min(plot_d$imageY)
    p1 <- ggplot2::ggplot(mapping = ggplot2::aes(minX:maxX, minX:maxY))
  }

  spe <- spe[features, ]
  plot_d <- SpatialExperiment::spatialCoords(spe) |>
    as.data.frame() |>
    dplyr::mutate(imageY = maxY - imageY)
  plot_d <- cbind(plot_d, as.matrix(t(SummarizedExperiment::assay(spe, assay_type))))
  colnames(plot_d) <- c('imageX', 'imageY', features)
  p1 + 
    scatterpie::geom_scatterpie(
      aes(x = imageX, y = imageY), 
      data = plot_d, cols = features, pie_scale = 0.6, color = NA) +
    coord_fixed() +
    theme_void()
}

#' Plot spatial pie chart of isoforms
#'
#' @description This function plots a spatial pie chart for given features.
#' @param spe The SpatialExperiment object.
#' @param isoforms The isoforms to plot.
#' @param assay_type The assay that contains the given features. E.g. 'counts', 'logcounts'.
#' @return A ggplot object.
#' @importFrom cowplot plot_grid
#' @export
plot_spatial_isoform <- function(spe, isoforms, assay_type = 'counts') {
  isoform_plot <- plot_isoforms(spe, transcript_ids = isoforms)
  pie_plot <- plot_spatial_pie(spe, isoforms, assay_type)
  cowplot::plot_grid(isoform_plot, pie_plot, ncol = 1)
}
