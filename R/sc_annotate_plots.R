#' Combine SCE
#'
#' Combine long- and short-read SingleCellExperiment objects
#'
#' @details Takes the long-read SCE object from the long-read pipeline and the short-read SCE object,
#' creates a \code{MultiAssayExperiment} object with the two \code{SingleCellExperiment} objects.
#' Cells with duplicated barcodes are removed from the larger library.
#'
#' @param short_read_large The SCE object, or path to the HDF5 file, or folder containing the matrix file,
#' corresponding to the larger short-read sample
#' @param short_read_small The SCE object, or path to the HDF5 file, or folder containing the matrix file,
#' corresponding to the smaller short-read sample
#' @param long_read_sce The SCE object of the transcript counts, from the long-read pipelines.
#' @param remove_duplicates determines whether cells with duplicated barcodes aer kept in the smaller library (
#' they are always removed from the larger library)
#'
#' @return A \code{MultiAssayExperiment} object, with 'gene_counts' and 'transcript_counts' experiments.
#'
#' @examples
#' library(SingleCellExperiment)
#' a <- SingleCellExperiment(assays = list(counts = matrix(rpois(100, 5), ncol = 10)))
#' b <- SingleCellExperiment(assays = list(counts = matrix(rpois(100, 5), ncol = 10)))
#' long_read <- SingleCellExperiment(assays = list(counts = matrix(rpois(100, 5), ncol = 10)))
#' colData(a)$Barcode <- paste0(1:10, '-1')
#' colData(b)$Barcode <- paste0(8:17, '-1')
#' colnames(long_read) <- as.character(2:11)
#' rownames(a) <- as.character(101:110)
#' rownames(b) <- as.character(103:112)
#' rownames(long_read) <- as.character(1001:1010)
#' combine_sce(short_read_large = a, short_read_small = b, long_read_sce = long_read)
#'
#' @importFrom BiocGenerics cbind colnames rownames
#' @importFrom SingleCellExperiment SingleCellExperiment altExp altExp<- counts counts<-
#' @importFrom SummarizedExperiment rowData colData rowData<- colData<- rowRanges rowRanges<-
#' @importFrom DropletUtils read10xCounts
#' @importFrom MultiAssayExperiment MultiAssayExperiment experiments<- experiments
#' @export
#' @md
combine_sce <- function(sce_with_lr, sce_without_lr) {
  if (is.character(sce_without_lr)) {
    sce_without_lr <- read10xCounts(sce_without_lr)
  }

  if (!is.null(colData(sce_with_lr)$Barcode) & !is.null(colData(sce_without_lr)$Barcode)) {
    duplicated_barcodes <- intersect(colData(sce_with_lr)$Barcode, colData(sce_without_lr)$Barcode)
    sce_without_lr <- sce_without_lr[, !colData(sce_without_lr)$Barcode %in% duplicated_barcodes]
  } else if (!is.null(colnames(sce_with_lr)) & !is.null(colnames(sce_without_lr))) {
    duplicated_barcodes <- intersect(colnames(sce_with_lr), colnames(sce_without_lr))
    sce_without_lr <- sce_without_lr[, !colnames(sce_without_lr) %in% duplicated_barcodes]
  }

  gene_intersect <- intersect(rownames(sce_with_lr), rownames(sce_without_lr))
  sce_with_lr <- sce_with_lr[gene_intersect,]
  sce_without_lr <- sce_without_lr[gene_intersect,]

  combined_sce <- SingleCellExperiment(assays = list(counts = cbind(counts(sce_with_lr), counts(sce_without_lr))))

  # altExp with NA values for sce_without_lr
  tr_counts <- cbind(
    counts(altExp(sce_with_lr, "transcript")),
    matrix(NA, nrow = nrow(altExp(sce_with_lr, "transcript")), ncol = ncol(sce_without_lr))
  )
  colnames(tr_counts) <- c(colnames(sce_with_lr), colnames(sce_without_lr))
  tr_sce <- SingleCellExperiment(assays = list(counts = tr_counts))

  rowData(tr_sce) <- rowData(altExp(sce_with_lr, "transcript"))
  rowRanges(tr_sce) <- rowRanges(altExp(sce_with_lr, "transcript"))
  rowData(combined_sce) <- rowData(sce_with_lr)
  rowRanges(combined_sce) <- rowRanges(sce_with_lr)

  altExp(combined_sce, "transcript") <- tr_sce
  combined_sce
}

sc_isoform_impute <- function(combined_sce, dimred = "PCA") {
  if (!"transcript" %in% names(altExps(combined_sce))) {
    stop("transcript counts not found")
  }
  if (!"logcounts" %in% names(assays(combined_sce))) {
    combined_sce <- scuttle::logNormCounts(combined_sce)
  }
  if (! "PCA" %in% reducedDimNames(combined_sce)) {
    combined_sce <- scater::runPCA(combined_sce)
  }
  if (!"UMAP" %in% reducedDimNames(combined_sce) & dimred == "UMAP") {
    combined_sce <- scater::runUMAP(combined_sce)
  }
  snn <- scran::buildSNNGraph(combined_sce, use.dimred = dimred)
  snn_mat <- igraph::as_adjacency_matrix(snn, attr = "weight")
  diag(snn_mat) <- ceiling(max(snn_mat))

  imputed_count <- sc_impute_expression(counts(altExp(combined_sce, "transcript")), snn_mat)
  assay(altExp(combined_sce, "transcript"), 'logcounts') <- imputed_count
  combined_sce
}

sc_impute_expression <- function(mtx, distance_matrix) {

  cat("Imputing transcript counts ...\n")
  expr_impute <- matrix(0, nrow(mtx), ncol(mtx))
  rownames(expr_impute) <- rownames(mtx)
  colnames(expr_impute) <- colnames(mtx)

  cols_with_counts <- apply(mtx, 2, function(x){!any(is.na(x))})
  expr_impute[, cols_with_counts] <- scuttle::normalizeCounts(mtx[, cols_with_counts])
  expr_impute <- expr_impute %*% distance_matrix
  colnames(expr_impute) <- colnames(mtx)

  expr_impute <- t(t(expr_impute)/colSums(distance_matrix))

  expr_impute
}


#' @importFrom stats quantile
reduce_quantile <- function(x, q = 0.05) {
  x[x < quantile(x, q, na.rm = TRUE)] <- quantile(x, q, na.rm = TRUE)
  x[x > quantile(x, 1 - q, na.rm = TRUE)] <- quantile(x, 1 - q, na.rm = TRUE)
  x
}


#' FLAMES UMAP plots
#'
#' Plot expression UMAPs of top n isoforms of a gene 
#'
#' @details
#' This function takes the combined \code{MultiAssayExperiment} object from 
#' \code{cexample("MultiAssayExperiment")ombine_sce} and plots UMAPs for each isoform of \code{gene}, where cells 
#' are colored by expression levels. When \code{grided = TRUE}, the UMAPs are combined into a grid, along with the isoforms' visualization along genomic coordinates. Produces a single UMAP with isoform expressions colored by \code{col_low} and \code{col_high} when \code{grided = FALSE}.
#'
#' @param gene The gene symbol of interest.
#' @param multiAssay The \code{MultiAssayExperiment} object from \code{combine_sce()}.
#' @param impute Whether to impute expression levels for cells without transcript counts
#' @param grided Wheter to produce multiple UMAP plots, with each showing expression level for an isoform, to allow plotting more than 2 isoforms.
#' @param n_isoforms The number of expressed isoforms to keep. \code{n_isoforms} > 2 requires \code{girded = TRUE}
#' @param transcript_ids specify the transcript ids instead of selecting the top \code{n_isoforms}
#' @param n_pcs The number of principal components to generate.
#' @param col_low Color for cells with low expression levels in UMAPs.
#' @param col_mid Color for cells with intermediate expression levels in UMAPs.
#' @param col_high Color for cells with high expression levels in UMAPs.
#'
#' @return a \code{ggplot} object of the UMAP(s)
#' 
#'
#' @examples
#' combined_sce <- combine_sce(
#'     short_read_large = scmixology_lib90,
#'     short_read_small = scmixology_lib10,
#'     long_read_sce = scmixology_lib10_transcripts,
#'     remove_duplicates = FALSE)
#'
#' sc_umap_expression(gene = "ENSG00000108107", multiAssay = combined_sce)
#' 
#' @importFrom MultiAssayExperiment experiments<- experiments
#' @importFrom cowplot get_legend plot_grid
#' @importFrom igraph as_adjacency_matrix
#' @importFrom scran buildSNNGraph fixedPCA getTopHVGs modelGeneVar
#' @importFrom stats median
#' @importFrom ggbio autoplot geom_alignment
#' @importFrom ggplot2 aes coord_polar element_blank element_line element_text 
#' @importFrom ggplot2 geom_bar geom_histogram geom_line geom_point geom_text 
#' @importFrom ggplot2 ggplot ggtitle labs margin position_stack 
#' @importFrom ggplot2 scale_colour_gradient2 theme theme_bw xlab ylab xlim
#' @export
#' @md
sc_umap_expression <- function(gene, multiAssay, impute = FALSE, grided = TRUE, n_isoforms = 4, transcript_ids,
  n_pcs = 40, col_low = "#313695", col_mid = "#FFFFBF", col_high = "#A50026") {
}

#' FLAMES heetmap plots
#'
#' Plot expression heatmap of top n isoforms of a gene 
#'
#' @details
#' This function takes the combined \code{MultiAssayExperiment} object from 
#' \code{combine_sce} and plots an expression heatmap with the isoform alignment visualisations.
#'
#' @param gene The gene symbol of interest.
#' @param multiAssay The \code{MultiAssayExperiment} object from \code{combine_sce()}.
#' @param impute Whether to impute expression levels for cells without transcript counts
#' @param n_isoforms The number of expressed isoforms to keep. 
#' @param transcript_ids specify the transcript ids instead of selecting the top \code{n_isoforms}
#' @param isoform_legend_width The width of isoform legends in heatmaps, in \code{cm}.
#' @param n_pcs The number of principal components to generate.
#' @param col_low Color for cells with low expression levels in UMAPs.
#' @param col_mid Color for cells with intermediate expression levels in UMAPs.
#' @param col_high Color for cells with high expression levels in UMAPs.
#' @param color_quantile The lower and upper expression quantile to be displayed bewteen \code{col_low} and \code{col_high}, e.g. with \code{color_quantile = 0.95}, cells with expressions higher than 95% of other cells will all be shown in \code{col_high}, and cells with expression lower than 95% of other cells will all be shown in \code{col_low}.
#'
#' @return a \code{ggplot} object of the heatmap
#' 
#' @examples
#' combined_sce <- combine_sce(
#'     short_read_large = scmixology_lib90,
#'     short_read_small = scmixology_lib10,
#'     long_read_sce = scmixology_lib10_transcripts,
#'     remove_duplicates = FALSE)
#' sc_heatmap_expression(gene = "ENSG00000108107", multiAssay = combined_sce)
#' 
#' @importFrom ComplexHeatmap AnnotationFunction Heatmap columnAnnotation rowAnnotation
#' @importFrom MultiAssayExperiment experiments<- experiments
#' @importFrom grid unit viewport
#' @importFrom gridExtra grid.arrange
#' @importFrom RColorBrewer brewer.pal
#' @importFrom circlize colorRamp2
#' @importFrom stats quantile
#' @importFrom igraph as_adjacency_matrix
#' @importFrom scran buildSNNGraph fixedPCA getTopHVGs modelGeneVar
#' @importFrom ggbio autoplot geom_alignment
#' @importFrom ggplot2 aes coord_polar element_blank element_line element_text 
#' @importFrom ggplot2 geom_bar geom_histogram geom_line geom_point geom_text 
#' @importFrom ggplot2 ggplot ggtitle labs margin position_stack 
#' @importFrom ggplot2 scale_colour_gradient2 theme theme_bw xlab ylab xlim
#' 
#' @export
#' @md
plot_isoform_heatmap <- function(sce, gene_id, transcript_ids, n = 4,
  isoform_legend_width = 7, col_low = "#313695", col_mid = "#FFFFBF", col_high = "#A50026", color_quantile = 0.95) {

  transcript_ids <- get_top_transcript_ids(sce, gene_id, transcript_ids, n)
  sce <- sce[match(transcript_ids, rowData(sce)$transcript_id),]
  legends_heatmap <- plot_top_isoforms(sce, gene_id, transcript_ids, n, format = 'list')

  group_annotation <- function(x) {
    n <- length(unique(x))
    if (n > 11 || n < 2) {
      column_anno <- HeatmapAnnotation(x)
    } else if (n == 2) {
      palette <- brewer.pal(n = 3, name = "PuRd")[c(1, 3)]
      names(palette) <- levels(x)
      column_anno <- HeatmapAnnotation(
        group = x,
        col = list(group = palette)
      )
    } else {
      palette <- brewer.pal(n, name = "Spectral")
      names(palette) <- unique(x)
      column_anno <- HeatmapAnnotation(
        group = x,
        col = list(group = palette)
      )
    }
    return(column_anno)
  }

  sce <- sce[, stats::hclust(stats::dist(t(logcounts(sce))))$order]
  expr_heatmap <- logcounts(sce) |>
    scale()

  expr_color_mapping <- function(expr_matrix) {
    if (color_quantile > 1 || color_quantile < 0) {
      color_quantile <- 0.95
    }
    return(colorRamp2(c(quantile(expr_matrix, 1 - color_quantile, na.rm = TRUE),
      0, quantile(expr_matrix, color_quantile, na.rm = TRUE)), c(col_low,
      col_mid, col_high)))
  }

  isoform_annotation <- AnnotationFunction(fun = function(index) {
    gridExtra::grid.arrange(grobs = legends_heatmap, ncol = 1, vp = viewport(),
      newpage = FALSE)
  }, which = "row", width = unit(isoform_legend_width, "cm"), n = length(legends_heatmap),
    subsettable = FALSE)

  return(
    Heatmap(expr_heatmap,
      name = "scaled expression",
      cluster_rows = FALSE, cluster_columns = FALSE, use_raster = FALSE,
      show_column_names = FALSE, show_row_names = FALSE, 
      # https://www.r-bloggers.com/2017/02/use-switch-instead-of-ifelse-to-return-a-null/
      top_annotation = switch(!is.null(colLabels(sce)), group_annotation(colLabels(sce))),
      left_annotation = rowAnnotation(isoform = isoform_annotation, annotation_name_rot = 0),
      col = expr_color_mapping(expr_heatmap)
    )
  )
}

get_top_transcript_ids <- function(sce, gene_id, transcript_ids, n){
  if (missing(transcript_ids)) {
    sce <- sce[rowData(sce)$gene_id == gene_id, ]
    sce <- sce[order(rowSums(counts(sce)), decreasing = TRUE), ]
    if (nrow(sce) > n) {
      sce <- sce[1:n, ]
    }
  } else {
    # sce <- sce[rowData(sce)$transcript_id %in% transcript_ids, ]
    # sce <- sce[order(rowSums(counts(sce)), decreasing = TRUE), ]
    return(transcript_ids)
  }
  return(rowData(sce)$transcript_id)
}

plot_top_isoforms <- function(sce, gene_id, transcript_ids, n = 4, format = 'plot_grid') {
  transcript_ids <- get_top_transcript_ids(sce, gene_id, transcript_ids, n)
  sce <- sce[match(transcript_ids, rowData(sce)$transcript_id),]
  transcripts <- rowRanges(sce)
  x_range <- c(
    min(min(start(transcripts))),
    max(max(end(transcripts)))
  )

  chr <- lapply(seqnames(transcripts), function(x){as.character(unique(x))}) |>
    unlist() |>
    unique()
  if (length(chr) > 1) {
    stop("transcripts span multiple chromosomes, unable to plot")
  }

  plot_list <- list()
  for (i in seq_along(transcript_ids)) {
    p <- ggplot(transcripts[i]) +
      ggbio::geom_alignment(
        label = TRUE, range.geom = "rect",
        gap.geom = "arrow", utr.geom = "rect"
      ) + 
      xlim(x_range) + 
      theme_void()
    if (i == length(transcript_ids)) {
      p <- p +
        theme(
          axis.line.x = element_line(),
          axis.ticks.x = element_line(),
          axis.text.x = element_text(),
          axis.title.x = element_text()
        )
        # labs(x = chr)
        # coord_cartesian(xlim = x_range)
    }
    plot_list <- append(plot_list, list(p@ggplot))
  }
  if (format == 'list') {
    return(plot_list)
  }
  return(cowplot::plot_grid(plotlist = plot_list, ncol = 1))
}

plot_isoform_umap <- function(sce, gene_id, transcript_ids, n = 4, reduced_dim_name = "UMAP", 
  expr_func = function(x){SingleCellExperiment::logcounts(x)}, 
  col_low = "#313695", col_mid = "#FFFFBF", col_high = "#A50026", format = 'plot_grid', ...) {
  transcript_ids <- get_top_transcript_ids(sce, gene_id, transcript_ids, n)
  sce <- sce[match(transcript_ids, rowData(sce)$transcript_id),]
  isoform_plot <- plot_top_isoforms(sce, gene_id, transcript_ids, n)
  df <- data.frame(
    x = reducedDim(sce, reduced_dim_name)[,1],
    y = reducedDim(sce, reduced_dim_name)[,2]
  )

  umaps <- list()
  for (i in seq_along(transcript_ids)) {
    df$expr <- expr_func(sce)[i,]
    p <- ggplot(df) +
      geom_point(aes(x = x, y = y, col = expr), alpha = 0.7, size = 0.2) +
      labs(x = "Dim1", y = "Dim2", col = "expression") +
      scale_colour_gradient2(low = col_low, mid = col_mid, high = col_high,
        na.value = "grey", midpoint = median(df$expr)) +
      theme_bw()
    umaps <- append(umaps, list(p))
  }

  plot_list <- c(list(isoform_plot), umaps)
  if (format == 'list') {
    return(umaps)
  }
  return(cowplot::plot_grid(plotlist = plot_list, ...))
}
