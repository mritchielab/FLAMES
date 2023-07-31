#' FLAMES Annotated Plottings
#'
#' (deprecated) Plot isoform exons alignments for a given gene, along with UMAP showing expression levels. Superseded by \code{sc_umap_expression} and \code{sc_heatmap_expression}.
#'
#'
#' @details
#' This function takes a combined \code{MultiAssayExperiment} object containing 'gene_counts' and 'transcript_counts' experiments 
#' to generate a combined UMAP, the expression levels of isoforms (using long read data)
#' are then overlayed on top of the UMAP. SNN inference based on gene counts were performed to impute isoform expression
#' for cells in the larger library. The \code{MultiAssayExperiment} object should have a boolean columen named 'Lib_small' in 
#' its columen data file to idicate which subsample the cells are in. The \code{MultiAssayExperiment} object can be created using
#' the \code{combine_sce} function.
#'
#' @param multiAssay The \code{MultiAssayExperiment} object from \code{combine_sce()}.
#' @param gene The gene symbol of interest.
#' @param n_isoforms The number of expressed isoforms to keep.
#' @param n_pcs The number of principal components to generate.
#' @param cluster_annotation Path to the cluster annotation CSV (required for heatmap, if \code{cluster_annotation.csv} is not in \code{path} and \code{multiAssay$cell_type} does not exist)
#' @param return_multiAssay Whether to return the processed \code{MultiAssayExperiment} object.
#' @param heatmap_annotation_colors Name of color palette to use for cell group annotation in heatmaps, see [RColorBrewer::brewer.pal()] \cr
#' available diverging palettes are:\cr
#' \code{BrBG PiYG PRGn PuOr RdBu RdGy RdYlBu RdYlGn Spectral}\cr
#' when there are more than 11 groups, this argument will be ignored and random palettes will be generated.
#' @param isoform_legend_width The width of isoform legends in heatmaps, in \code{cm}.
#' @param col_low Color for cells with low expression levels in UMAPs.
#' @param col_mid Color for cells with intermediate expression levels in UMAPs.
#' @param col_high Color for cells with high expression levels in UMAPs.
#' @param heatmap_color_quantile Float; Expression levels higher than this quantile will all be shown with \code{col_high}.
#' Expression levels lower than 1 - \code{heatmap_color_quantile} will all be shown with \code{col_low};
#'
#' @return a list containing the combined UMAP, the isoform exon alignments and the UMAP with isoform expression levels.
#' @seealso
#' [combine_sce()] for combining gene count \code{SingleCellExperiment} object with transcript counts.
#'
#' @importFrom dplyr group_by summarise_at slice_max filter
#' @importFrom tidyr gather pivot_wider
#' @importFrom magrittr '%>%'
#' @importFrom S4Vectors DataFrame
#' @importFrom SingleCellExperiment SingleCellExperiment reducedDimNames logcounts altExp altExp<-
#' @importFrom SummarizedExperiment rowData colData rowData<- colData<-
#' @importFrom scuttle addPerCellQC addPerFeatureQC isOutlier logNormCounts
#' @importFrom scran getTopHVGs fixedPCA buildSNNGraph modelGeneVar
#' @importFrom scater runUMAP
#' @importFrom ggplot2 ggplot geom_point aes labs element_blank element_line element_text theme_bw theme scale_colour_gradient2 margin
#' @importFrom ggbio autoplot geom_alignment xlim
#' @importFrom GenomicRanges GRanges GRangesList
#' @importFrom rtracklayer import.gff3
#' @importFrom BiocGenerics cbind colnames rownames start end
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation columnAnnotation AnnotationFunction rowAnnotation
#' @importFrom igraph as_adjacency_matrix
#' @importFrom rtracklayer import
#' @importFrom Matrix t colSums
#' @importFrom cowplot plot_grid get_legend
#' @importFrom RColorBrewer brewer.pal
#' @importFrom circlize colorRamp2
#' @importFrom grid unit viewport
#' @importFrom gridExtra grid.arrange
#' @importFrom stats quantile median
#' @importFrom MultiAssayExperiment experiments experiments<-
#' @export
#' @md
sc_annotate_plots <- function(gene, multiAssay, cluster_annotation, n_isoforms = 4, n_pcs = 40, return_multiAssay = TRUE,
                              heatmap_annotation_colors = "BrBG", isoform_legend_width = 7, heatmap_color_quantile = 0.95, col_low = "#313695", col_mid = "#FFFFBF", col_high = "#A50026") {

  # Check args
  if (is.null(colData(multiAssay)$Lib_small)) {
    stop("Please provide colData(multiAssay)$Lib_small
          to specify the library of each cell
          Values should be either TRUE or FALSE.\n")
  }

  outputs <- list()

  ### Process SCEs
  flames_outdir <- experiments(multiAssay)$transcript_counts@metadata$OutputFiles$outdir
  experiments(multiAssay)$transcript_counts <- logNormCounts(experiments(multiAssay)$transcript_counts)
  if (!"mean" %in% names(rowData(experiments(multiAssay)$transcript_counts))) {
    experiments(multiAssay)$transcript_counts <- addPerFeatureQC(experiments(multiAssay)$transcript_counts) # needed for sorting
  }

  if (!("PCA" %in% reducedDimNames(experiments(multiAssay)$gene_counts)) || dim(experiments(multiAssay)$gene_counts@int_colData$reducedDims$PCA)[2] < n_pcs) {
    cat("Running PCA for experiments(multiAssay)$gene_counts ...\n")
    experiments(multiAssay)$gene_counts <- logNormCounts(experiments(multiAssay)$gene_counts)
    hvgs <- getTopHVGs(modelGeneVar(experiments(multiAssay)$gene_counts), n = 2000)
    experiments(multiAssay)$gene_counts <- fixedPCA(experiments(multiAssay)$gene_counts, rank = n_pcs, subset.row = hvgs)
    snn <- buildSNNGraph(experiments(multiAssay)$gene_counts, use.dimred = "PCA")
  } else if (dim(experiments(multiAssay)$gene_counts@int_colData$reducedDims$PCA)[2] > n_pcs) {
    message(paste0(c(
      "experiments(multiAssay)$gene_counts have ", dim(experiments(multiAssay)$gene_counts@int_colData$reducedDims$PCA)[2], " PCs, using the first ",
      n_pcs, " PCs since n_pcs set to ", n_pcs, "\n"
    )))
    experiments(multiAssay)$gene_counts@int_colData$reducedDims$tmp <- experiments(multiAssay)$gene_counts@int_colData$reducedDims$PCA[, 1:n_pcs]
    snn <- buildSNNGraph(experiments(multiAssay)$gene_counts, use.dimred = "tmp")
    experiments(multiAssay)$gene_counts@int_colData$reducedDims$tmp <- NULL
  } else {
    snn <- buildSNNGraph(experiments(multiAssay)$gene_counts, use.dimred = "PCA")
  }

  # Distance matrix for imputation
  snn_mat <- as_adjacency_matrix(snn, attr = "weight")
  diag(snn_mat) <- ceiling(max(snn_mat))

  if (!("UMAP" %in% reducedDimNames(experiments(multiAssay)$gene_counts))) {
    cat("Running UMAP for experiments(multiAssay)$gene_counts ...\n")
    experiments(multiAssay)$gene_counts <- runUMAP(experiments(multiAssay)$gene_counts, dimred = "PCA")
  } else {
    cat("Skipping runUMAP...\n")
  }

  sample_lib <- factor(multiAssay[,,"gene_counts"]$Lib_small)
  levels(sample_lib) <- gsub("FALSE", "large", gsub("TRUE", "small", levels(sample_lib)))
  names(sample_lib) <- colnames(multiAssay[,,"gene_counts"])[["gene_counts"]]

  plot_umap <- ggplot() +
    geom_point(aes(
      x = experiments(multiAssay)$gene_counts@int_colData$reducedDims$UMAP[, 1],
      y = experiments(multiAssay)$gene_counts@int_colData$reducedDims$UMAP[, 2],
      col = sample_lib
    ),
    size = 0.02
    ) +
    labs(x = "Dim1", y = "Dim2", title = "Lib_all UMAP", color = "library")

  # row_meta columns: "transcript_id", "gene_id", "FSM_match", "mean", "detected", "gene_name"
  # rows: transcript_id
  row_meta <- as.data.frame(rowData(experiments(multiAssay)$transcript_counts))
  row_meta <- row_meta %>%
    dplyr::filter(gene_id == gene) %>%
    slice_max(n = n_isoforms, order_by = mean) # keep top n isoforms per gene
  if (dim(row_meta)[1] < 2) {
    stop(paste0("Alternative isoform not found for the given gene: ", gene))
  }

  # apply top n to sce
  tr_sce_multi <- experiments(multiAssay[rownames(row_meta),
        names(sample_lib)[which(sample_lib == "small")],
        "transcript_counts"])$transcript_counts 
  if (n_isoforms > nrow(row_meta)) {
    n_isoforms <- nrow(row_meta)
    message(paste0(c("Only ", n_isoforms, " isoforms found for this gene.\n")))
  }

  # Alignment info
  isoform_sel <- rowRanges(experiments(multiAssay)$transcript_counts[row_meta$FSM_match, ])
  # isoform_sel <- isoform_sel[row_meta$FSM_match] # Order by expression levels
  if (length(isoform_sel) == 2) {
    fill_by_isoform <- c(rep(col_low, length(isoform_sel[[1]])), rep(col_high, length(isoform_sel[[2]])))
    plot_isoforms <- ggbio::autoplot(isoform_sel, label = TRUE, fill = fill_by_isoform) +
      theme_bw() + theme(
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")
      )
  } else {
    plot_isoforms <- ggbio::autoplot(isoform_sel, label = TRUE) +
      theme_bw() + theme(
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")
      )
    # ggbio::geom_alignment does not follow the order in isoform_sel
  }

  # Could try using stat_stepping() to get the order of isoforms
  # https://github.com/lawremi/ggbio/issues/149
  tr_order <- rev(names(isoform_sel)[match(1:n_isoforms, plot_isoforms$layers[[length(plot_isoforms$layers)]]$data$stepping)])
  # isoform_sel <- isoform_sel[tr_order]

  legends_heatmap <- list()
  for (i in names(isoform_sel)) {
    p <- ggplot(isoform_sel[[i]]) +
      geom_alignment(
        label = FALSE, range.geom = "rect",
        gap.geom = "arrow", utr.geom = "rect"
      ) +
      xlim(c(min(min(start(isoform_sel))), max(max(end(isoform_sel))))) +
      theme(
        axis.line = element_blank(), axis.text.x = element_blank(),
        axis.text.y = element_blank(), axis.ticks = element_blank(),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        legend.position = "none", panel.background = element_blank(),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), plot.background = element_blank(),
        plot.margin = margin(20, 0, 20, 0), plot.title = element_text(hjust = 0.5)
      ) +
      labs(title = i)
    legends_heatmap[[i]] <- p@ggplot
  }
  legends_heatmap <- legends_heatmap[tr_order]

  ### impute transcript counts
  cat("Imputing transcript counts ...\n")
  expr <- logcounts(tr_sce_multi)
  expr_all <- matrix(0, nrow(expr), length(sample_lib)) # dim: transcripts x cells, initialised to 0
  rownames(expr_all) <- rownames(expr)
  colnames(expr_all) <- names(sample_lib)

  # expr <- expr[, colnames(expr) %in% colnames(multiAssay)]
  expr_all[, colnames(expr)] <- expr
  expr_all <- expr_all %*% snn_mat
  colnames(expr_all) <- names(sample_lib)

  expr_all <- t(t(expr_all) / colSums(snn_mat))
  expr_all <- scale(expr_all)

  reduce_quantile <- function(x, q = 0.05) {
    x[x < quantile(x, q, na.rm = TRUE)] <- quantile(x, q, na.rm = TRUE)
    x[x > quantile(x, 1 - q, na.rm = TRUE)] <- quantile(x, 1 - q, na.rm = TRUE)
    x
  }
  expr_all <- t(apply(expr_all, 1, reduce_quantile))

  # Expressions without imputation (Lib_20)
  expr_20 <- matrix(0, nrow(expr), sum(sample_lib == "small"))
  rownames(expr_20) <- rownames(expr)
  colnames(expr_20) <- names(sample_lib[sample_lib == "small"])
  expr_20[, colnames(expr)] <- expr
  expr_20 <- scale(expr_20)
  expr_20 <- t(apply(expr_20, 1, reduce_quantile))

  # scale() will produce NaN values for cell with 0 counts for all isoforms
  # using na.value = "grey" for UMAPs

  ### Plot with UMAP from lib_all
  cat("Plotting expression UMAPs ...\n")
  umap_all <- as.data.frame(experiments(multiAssay)$gene_counts@int_colData$reducedDims$UMAP)
  colnames(umap_all) <- c("x", "y")

  umap_80 <- umap_all[sample_lib == "large", ]
  colnames(umap_80) <- c("x", "y")

  umap_20 <- umap_all[sample_lib == "small", ]
  colnames(umap_20) <- c("x", "y")

  if (n_isoforms == 2) {
    # issue: Rps24.png
    umap_20$expr <- "grey"
    umap_all$expr <- "grey"
    umap_20$expr[expr_20[names(isoform_sel)[1], ] > 0] <- col_low
    umap_all$expr[expr_all[names(isoform_sel)[1], ] > 0] <- col_low
    umap_20$expr[expr_20[names(isoform_sel)[1], ] < 0] <- col_high
    umap_all$expr[expr_all[names(isoform_sel)[1], ] < 0] <- col_high
    plot_expression_umaps <- ggplot() +
      geom_point(data = umap_80, aes(x = x, y = y), alpha = 0.2, size = 0.2, col = "grey", show.legend = FALSE) +
      geom_point(data = umap_20, aes(x = x, y = y), col = umap_20$expr, alpha = 0.7, size = 0.7) +
      labs(x = "Dim1", y = "Dim2", col = "dominant isoform") +
      theme_bw() +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.text = element_blank()
      )
    plot_expression_umaps_impute <- ggplot() +
      geom_point(data = umap_all, aes(x = x, y = y), col = umap_all$expr, alpha = 0.7, size = 0.2) +
      labs(x = "Dim1", y = "Dim2", col = "dominant isoform") +
      theme_bw() +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.text = element_blank()
      )
    combined_isoform_plot <- plot_grid(plot_isoforms@ggplot, plot_expression_umaps, ncol = 1, rel_heights = c(1, 2.5))
    combined_isoform_plot_impute <- plot_grid(plot_isoforms@ggplot, plot_expression_umaps_impute, ncol = 1, rel_heights = c(1, 2.5))
  } else {
    umap_20 <- cbind(t(expr_20), umap_20)
    umap_all <- cbind(t(expr_all), umap_all)
    # ggplot2 lazy evaluation with for loops
    #   may need to unquote idx with !! operator
    plot_idx <- function(idx) {
      p <- ggplot() +
        geom_point(data = umap_80, aes(x = x, y = y), alpha = 0.2, size = 0.2, col = "grey", show.legend = FALSE) +
        geom_point(data = umap_20, aes(x = x, y = y, col = umap_20[, idx]), alpha = 0.7, size = 0.7) +
        labs(x = "Dim1", y = "Dim2", col = "scaled expression", title = colnames(umap_20)[idx]) +
        scale_colour_gradient2(low = col_low, mid = col_mid, high = col_high, na.value = "grey", midpoint = median(umap_20[, idx], na.rm = TRUE)) +
        theme_bw() +
        theme(
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = "none"
        )
      p
    }
    plot_idx_impute <- function(idx) {
      p <- ggplot() +
        geom_point(data = umap_all, aes(x = x, y = y, col = umap_all[, idx]), alpha = 0.7, size = 0.2) +
        labs(x = "Dim1", y = "Dim2", col = "scaled expression", title = colnames(umap_all)[idx]) +
        scale_colour_gradient2(low = col_low, mid = col_mid, high = col_high, na.value = "grey", midpoint = median(umap_all[, idx], na.rm = TRUE)) +
        theme_bw() +
        theme(
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = "none"
        )
      p
    }
    plot_expression_umaps <- lapply(1:n_isoforms, plot_idx)
    plot_expression_umaps_impute <- lapply(1:n_isoforms, plot_idx_impute)
    plot_expression_umaps <- plot_expression_umaps[match(tr_order, colnames(umap_20)[1:n_isoforms])]
    plot_expression_umaps_impute <- plot_expression_umaps_impute[match(tr_order, colnames(umap_all)[1:n_isoforms])]


    legend <- ggplot(data.frame(x = 1:2, y = 1:2)) +
      geom_point(aes(x = x, y = y, col = c(-1, 1))) +
      labs(col = "scaled expression") +
      scale_colour_gradient2(
        low = col_low, mid = col_mid, high = col_high,
        breaks = c(-0.8, 0.8),
        labels = c("Low", "High")
      )
    legend <- get_legend(legend + theme(legend.box.margin = margin(0, 0, 0, 12)))
    combined_isoform_plot <- plot_grid(plotlist = c(list(annotation = plot_isoforms@ggplot), plot_expression_umaps))
    combined_isoform_plot <- plot_grid(combined_isoform_plot, legend, rel_widths = c(3, .4))
    combined_isoform_plot_impute <- plot_grid(plotlist = c(list(annotation = plot_isoforms@ggplot), plot_expression_umaps_impute))
    combined_isoform_plot_impute <- plot_grid(combined_isoform_plot_impute, legend, rel_widths = c(3, .4))
  }

  #  Cluster info
  if ((!is.null(flames_outdir) && file.exists(file.path(flames_outdir, "cluster_annotation.csv"))) ||
    (!missing(cluster_annotation) && file.exists(cluster_annotation)) ||
    !is.null(multiAssay$cell_type)) {
    cat("Plotting heatmaps ...\n")
    if (file.exists(file.path(flames_outdir, "cluster_annotation.csv"))) {
      cluster_barcode <- read.csv(file.path(flames_outdir, "cluster_annotation.csv"), stringsAsFactors = FALSE)
      cluster_barcode <- cluster_barcode[, c("barcode_seq", "groups")]
    } else if (!missing(cluster_annotation) && file.exists(cluster_annotation)) {
      cluster_barcode <- read.csv(cluster_annotation, stringsAsFactors = FALSE)
      cluster_barcode <- cluster_barcode[, c("barcode_seq", "groups")]
    } else {
      cluster_barcode <- data.frame(barcode_seq = names(sample_lib), groups = multiAssay[,names(sample_lib),]$cell_type)
    }
    rownames(cluster_barcode) <- cluster_barcode[, "barcode_seq"]

    outputs[["umap_clusters"]] <- ggplot(umap_all) +
      geom_point(aes(
        x = x, y = y,
        col = factor(cluster_barcode[names(sample_lib), "groups"])
      ),
      size = 0.2, alpha = 0.7
      ) +
      labs(x = "dim1", y = "dim2", title = "lib_all umap", color = "cluster")

    # Heatmaps

    # Error in hclust(x) : NA/NaN/Inf in foreign function call (arg 10)
    # removing NaN cells from expr_20 and expr_all
    expr_20 <- expr_20[, !apply(expr_20, 2, function(x) {
      any(is.na(x))
    })]
    expr_all <- expr_all[, !apply(expr_all, 2, function(x) {
      any(is.na(x))
    })]

    cell_order_20 <- stats::hclust(stats::dist(t(expr_20)))$order
    cell_order_all <- stats::hclust(stats::dist(t(expr_all)))$order
    cluster_barcode_20 <- cluster_barcode[colnames(expr_20)[cell_order_20], "groups", drop = FALSE]
    cluster_barcode_all <- cluster_barcode[colnames(expr_all)[cell_order_all], "groups", drop = FALSE]


    group_annotation <- function(cluster_barcode) {
      n <- length(table(cluster_barcode))
      if (n > 11 || n < 2) {
        column_anno <- columnAnnotation(df = cluster_barcode)
      } else if (n == 2) {
        palette <- brewer.pal(n = 3, name = heatmap_annotation_colors)[c(1, 3)]
        names(palette) <- names(table(cluster_barcode))
        column_anno <- columnAnnotation(
          df = cluster_barcode,
          col = list(groups = palette)
        )
      } else {
        palette <- brewer.pal(n, name = heatmap_annotation_colors)
        names(palette) <- names(table(cluster_barcode))
        column_anno <- columnAnnotation(
          df = cluster_barcode,
          col = list(groups = palette)
        )
      }
      return(column_anno)
    }

    expr_color_mapping <- function(expr_matrix) {
      if (heatmap_color_quantile > 1 || heatmap_color_quantile < 0) {
        heatmap_color_quantile <- 0.95
      }
      return(colorRamp2(
        c(quantile(expr_matrix, 1 - heatmap_color_quantile, na.rm = TRUE), 0, quantile(expr_matrix, heatmap_color_quantile, na.rm = TRUE)),
        c(col_low, col_mid, col_high)
      ))
    }

    isoform_annotation <- AnnotationFunction(
      fun = function(index) {
        gridExtra::grid.arrange(grobs = legends_heatmap, ncol = 1, vp = viewport(), newpage = FALSE)
      },
      which = "row",
      width = unit(isoform_legend_width, "cm"),
      n = length(legends_heatmap),
      subsettable = FALSE
    )

    expr_20 <- expr_20[tr_order, cell_order_20]
    expr_all <- expr_all[tr_order, cell_order_all]


    outputs[["heatmap"]] <- Heatmap(expr_20,
      name = "scaled expression",
      cluster_rows = FALSE, cluster_columns = FALSE, use_raster = FALSE, show_column_names = FALSE,
      top_annotation = group_annotation(cluster_barcode_20), col = expr_color_mapping(expr_20)
    )

    outputs[["heatmap_impute"]] <- Heatmap(expr_all,
      name = "scaled expression",
      cluster_rows = FALSE, cluster_columns = FALSE, use_raster = FALSE, show_column_names = FALSE,
      top_annotation = group_annotation(cluster_barcode_all), col = expr_color_mapping(expr_all)
    )

    outputs[["combined_heatmap"]] <- Heatmap(expr_20,
      name = "scaled expression",
      cluster_rows = FALSE, cluster_columns = FALSE, use_raster = FALSE,
      show_column_names = FALSE, show_row_names = FALSE, top_annotation = group_annotation(cluster_barcode_20),
      left_annotation = rowAnnotation(isoform = isoform_annotation, annotation_name_rot = 0),
      col = expr_color_mapping(expr_20)
    )

    outputs[["combined_heatmap_impute"]] <- Heatmap(expr_all,
      name = "scaled expression",
      cluster_rows = FALSE, cluster_columns = FALSE, use_raster = FALSE,
      show_column_names = FALSE, show_row_names = FALSE, top_annotation = group_annotation(cluster_barcode_all),
      left_annotation = rowAnnotation(isoform = isoform_annotation, annotation_name_rot = 0),
      col = expr_color_mapping(expr_all)
    )
  } else {
    cat("Cluster annotation not found, heatmaps skipped.\n")
  }

  outputs[["plot_umap"]] <- plot_umap
  outputs[["plot_isoforms"]] <- plot_isoforms
  outputs[["plot_expression_umaps"]] <- plot_expression_umaps
  outputs[["plot_expression_umaps_impute"]] <- plot_expression_umaps_impute
  outputs[["combined_isoform_plot"]] <- combined_isoform_plot
  outputs[["combined_isoform_plot_impute"]] <- combined_isoform_plot_impute


  if (return_multiAssay) {
    outputs[["multiAssay"]] <- multiAssay
  }

  return(outputs)
}

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
#' @importFrom MultiAssayExperiment MultiAssayExperiment
#' @export
#' @md
combine_sce <- function(short_read_large, short_read_small, long_read_sce, remove_duplicates = FALSE) {
  if (is.character(short_read_large)) {
    short_read_large <- read10xCounts(short_read_large)
  }
  stopifnot(is(short_read_large, "SingleCellExperiment"))
  if (is.character(short_read_small)) {
    short_read_small <- read10xCounts(short_read_small)
  }
  stopifnot(is(short_read_small, "SingleCellExperiment"))

  if (is.null(colnames(short_read_large)) || is.null(colnames(short_read_small))) {
    if ("Barcode" %in% names(colData(short_read_large)) && "Barcode" %in% names(colData(short_read_small))) {
      colnames(short_read_large) <- gsub("-1$", "", colData(short_read_large)$Barcode)
      colnames(short_read_small) <- gsub("-1$", "", colData(short_read_small)$Barcode)
    } else {
      stop("Error getting cell barcode.")
    }
  }

  # Always remove cells with duplicated barcodes from the larger library
  # remove_duplicates determines whether they aer kept in the smaller library
  dup_bc <- colnames(short_read_small[, colnames(short_read_small) %in% colnames(short_read_large)])
  short_read_large <- short_read_large[, !(colnames(short_read_large) %in% dup_bc)]
  if (remove_duplicates) {
    short_read_small <- short_read_small[, !(colnames(short_read_small) %in%
      dup_bc)]
  }

  gene_intersection <- intersect(rownames(short_read_large), rownames(short_read_small))
  short_read_all <- cbind(short_read_large[gene_intersection, ], short_read_small[gene_intersection,
    ])

  sampleMap <- rbind(DataFrame(assay = "gene_counts", primary = colnames(short_read_all),
    colname = colnames(short_read_all)), DataFrame(assay = "transcript_counts",
    primary = colnames(long_read_sce), colname = colnames(long_read_sce)))

  colLib <- DataFrame(row.names = unique(sampleMap$colname))
  colLib$Lib_small <- rownames(colLib) %in% colnames(short_read_small)

  return(MultiAssayExperiment::MultiAssayExperiment(list(gene_counts = short_read_all,
    transcript_counts = long_read_sce), sampleMap = sampleMap, colData = colLib))
}

sc_reduce_dims <- function(multiAssay, n_pcs = 40, n_hvgs = 2000) {

  experiments(multiAssay)$transcript_counts <- logNormCounts(experiments(multiAssay)$transcript_counts)
  if (!"mean" %in% names(rowData(experiments(multiAssay)$transcript_counts))) {
    experiments(multiAssay)$transcript_counts <- addPerFeatureQC(experiments(multiAssay)$transcript_counts)  # needed for sorting
  }

  if (!("PCA" %in% reducedDimNames(experiments(multiAssay)$gene_counts)) || dim(experiments(multiAssay)$gene_counts@int_colData$reducedDims$PCA)[2] <
    n_pcs) {
    cat("Running PCA for experiments(multiAssay)$gene_counts ...\n")
    experiments(multiAssay)$gene_counts <- logNormCounts(experiments(multiAssay)$gene_counts)
    hvgs <- getTopHVGs(modelGeneVar(experiments(multiAssay)$gene_counts), n = n_hvgs)
    experiments(multiAssay)$gene_counts <- fixedPCA(experiments(multiAssay)$gene_counts,
      rank = n_pcs, subset.row = hvgs)
  } else if (dim(experiments(multiAssay)$gene_counts@int_colData$reducedDims$PCA)[2] >
    n_pcs) {
    message(paste0(c("experiments(multiAssay)$gene_counts have ", dim(experiments(multiAssay)$gene_counts@int_colData$reducedDims$PCA)[2],
      " PCs, using the first ", n_pcs, " PCs since n_pcs set to ", n_pcs, "\n")))
    experiments(multiAssay)$gene_counts@int_colData$reducedDims$tmp <- experiments(multiAssay)$gene_counts@int_colData$reducedDims$PCA[,
      1:n_pcs]
    experiments(multiAssay)$gene_counts@int_colData$reducedDims$tmp <- NULL
  }

  if (!("UMAP" %in% reducedDimNames(experiments(multiAssay)$gene_counts))) {
    cat("Running UMAP for experiments(multiAssay)$gene_counts ...\n")
    experiments(multiAssay)$gene_counts <- runUMAP(experiments(multiAssay)$gene_counts,
      dimred = "PCA")
  } else {
    cat("Skipping runUMAP...\n")
  }

  return(multiAssay)
}

reduce_quantile <- function(x, q = 0.05) {
  x[x < quantile(x, q, na.rm = TRUE)] <- quantile(x, q, na.rm = TRUE)
  x[x > quantile(x, 1 - q, na.rm = TRUE)] <- quantile(x, 1 - q, na.rm = TRUE)
  x
}

sc_impute_expression <- function(expression_subsample, distance_matrix, sample_lib, scale_exprsion = TRUE) {

  cat("Imputing transcript counts ...\n")
  expr_all <- matrix(0, nrow(expression_subsample), length(sample_lib))  # dim: transcripts x cells, initialised to 0
  rownames(expr_all) <- rownames(expression_subsample)
  colnames(expr_all) <- names(sample_lib)

  expr_all[, colnames(expression_subsample)] <- expression_subsample
  expr_all <- expr_all %*% distance_matrix
  colnames(expr_all) <- names(sample_lib)

  expr_all <- t(t(expr_all)/colSums(distance_matrix))

  if (scale_exprsion) {
  expr_all <- scale(expr_all)
  expr_all <- t(apply(expr_all, 1, reduce_quantile))
  }
  expr_all
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
#' source(system.file("examples/scmixology1.R", package = "FLAMES"))
#' scm_lib80 <- scuttle::addPerCellQC(scm_lib80)
#' scm_lib20 <- scuttle::addPerCellQC(scm_lib20)
#' qc_80 <- scuttle::quickPerCellQC(colData(scm_lib80))
#' qc_20 <- scuttle::quickPerCellQC(colData(scm_lib20))
#' qc_80$discard[scm_lib80$sum < 10000] <- TRUE
#' qc_20$discard[scm_lib20$sum < 20000] <- TRUE
#' 
#' combined_sce <- combine_sce(
#'     short_read_large = scm_lib80[,!qc_80$discard],
#'     short_read_small = scm_lib20[,!qc_20$discard],
#'     long_read_sce = scm_lib20_transcripts,
#'     remove_duplicates = FALSE)
#'
#' sc_umap_expression(gene = "ENSG00000108107", multiAssay = combined_sce)
#' 
#' @export
#' @md
sc_umap_expression <- function(gene, multiAssay, impute = FALSE, grided = TRUE, n_isoforms = 4, transcript_ids,
  n_pcs = 40, col_low = "#313695", col_mid = "#FFFFBF", col_high = "#A50026") {

  multiAssay <- sc_reduce_dims(multiAssay, n_pcs)

  sample_lib <- factor(multiAssay[, , "gene_counts"]$Lib_small)
  levels(sample_lib) <- gsub("FALSE", "large", gsub("TRUE", "small", levels(sample_lib)))
  names(sample_lib) <- colnames(multiAssay[, , "gene_counts"])[["gene_counts"]]

  if (!grided) {
    if (!missing(transcript_ids) && length(transcript_ids) > 2) {
      stop("grided = FALSE can only plot 2 isoforms, consider setting it to TRUE\n")
    } 
    if (!missing(n_isoforms) && n_isoforms > 2) {
      stop("grided = FALSE can only plot 2 isoforms, consider setting n_isoforms = 2 or
        grided = TRUE\n")} 
  }

  # row_meta columns: 'transcript_id', 'gene_id', 'FSM_match', 'mean',
  # 'detected', 'gene_name' rows: transcript_id
  row_meta <- as.data.frame(rowData(experiments(multiAssay)$transcript_counts))
  if (!missing(transcript_ids) && !is.null(transcript_ids)) {
    if (!all(transcript_ids %in% rownames(row_meta))) {
      message(paste0(c("transcript_id(s) ", paste(transcript_ids[!(transcript_ids %in%
        rownames(transcript_id))], sep = " "), " not found\n")))
    }
    row_meta <- row_meta %>%
      dplyr::filter(transcript_id %in% transcript_ids)
    stopifnot(`no transcript_id match` = dim(row_meta)[1] == 0)
  } else {
    row_meta <- row_meta %>%
      dplyr::filter(gene_id == gene) %>%
      slice_max(n = n_isoforms, order_by = mean)  # keep top n isoforms per gene
    if (dim(row_meta)[1] < 2) {
      stop(paste0("Alternative isoform not found for the given gene: ", gene))
    }
  }

  # apply selection to sce
  tr_sce_multi <- experiments(multiAssay[rownames(row_meta), names(sample_lib)[which(sample_lib ==
    "small")], "transcript_counts"])$transcript_counts
  if (n_isoforms > nrow(row_meta)) {
    n_isoforms <- nrow(row_meta)
    message(paste0(c("Only ", n_isoforms, " isoforms found for this gene.\n")))
  }

  # Isoform alignment visualization
  isoform_sel <- rowRanges(experiments(multiAssay)$transcript_counts[row_meta$FSM_match,
    ])
  # isoform_sel <- isoform_sel[row_meta$FSM_match] # Order by expression levels
  if (grided) {
    plot_isoforms <- ggbio::autoplot(isoform_sel, label = TRUE) + theme_bw() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
    # ggbio::geom_alignment does not follow the order in isoform_sel
  } else {
    fill_by_isoform <- c(rep(col_low, length(isoform_sel[[1]])), rep(col_high,
      length(isoform_sel[[2]])))
    plot_isoforms <- ggbio::autoplot(isoform_sel, label = TRUE, fill = fill_by_isoform) +
      theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.background = element_blank(), axis.line = element_line(colour = "black"))
  }

  # Could try using stat_stepping() to get the order of isoforms
  # https://github.com/lawremi/ggbio/issues/149
  tr_order <- rev(names(isoform_sel)[match(1:n_isoforms, plot_isoforms$layers[[length(plot_isoforms$layers)]]$data$stepping)])
  # isoform_sel <- isoform_sel[tr_order]

  # expression levels
  if (impute) {
    snn <- buildSNNGraph(experiments(multiAssay)$gene_counts, use.dimred = "PCA")
    snn_mat <- as_adjacency_matrix(snn, attr = "weight")
    diag(snn_mat) <- ceiling(max(snn_mat))
    expr_all <- sc_impute_expression(logcounts(tr_sce_multi), snn_mat, sample_lib)
  }
  # Expressions without imputation (Lib_20)
  expr_20 <- matrix(0, dim(tr_sce_multi)[1], sum(sample_lib == "small"))
  rownames(expr_20) <- rownames(tr_sce_multi)
  colnames(expr_20) <- names(sample_lib[sample_lib == "small"])
  expr_20[, colnames(tr_sce_multi)] <- logcounts(tr_sce_multi)
  expr_20 <- scale(expr_20)
  expr_20 <- t(apply(expr_20, 1, reduce_quantile))

  # scale() will produce NaN values for cell with 0 counts for all isoforms
  # using na.value = 'grey' for UMAPs

  ### Plot with UMAP from lib_all
  cat("Plotting expression UMAPs ...\n")
  umap_all <- as.data.frame(experiments(multiAssay)$gene_counts@int_colData$reducedDims$UMAP)
  colnames(umap_all) <- c("x", "y")

  umap_80 <- umap_all[sample_lib == "large", ]
  colnames(umap_80) <- c("x", "y")

  umap_20 <- umap_all[sample_lib == "small", ]
  colnames(umap_20) <- c("x", "y")

  umap_20 <- cbind(t(expr_20), umap_20)

  ggplot_theme_bw <- theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), panel.border = element_blank(), axis.line = element_line(colour = "black"),
    axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(),
    axis.ticks.y = element_blank(), legend.position = "none")

  if (!grided) {
    if (impute) {
      gradient_col <- expr_all[names(isoform_sel)[1], ]/(expr_all[names(isoform_sel)[1],
        ] + expr_all[names(isoform_sel)[2], ])
      plot_expression_umaps <- ggplot() + 
      geom_point(data = umap_all, aes(x = x, y = y, col = gradient_col), alpha = 0.7, size = 0.7) + 
      labs(x = "Dim1", y = "Dim2", col = "scaled imputed expression") +
      scale_colour_gradient2(low = col_low, mid = col_mid, high = col_high, na.value = "grey", midpoint = 0.5) + 
      ggplot_theme_bw
    } else {
      gradient_col <- expr_20[names(isoform_sel)[1], ]/(expr_20[names(isoform_sel)[1],
        ] + expr_20[names(isoform_sel)[2], ])
      plot_expression_umaps <- ggplot() + 
      geom_point(data = umap_80, aes(x = x, y = y), alpha = 0.2, size = 0.2, col = "grey", show.legend = FALSE) +
      geom_point(data = umap_20, aes(x = x, y = y, col = gradient_col), alpha = 0.7, size = 0.7) + 
      labs(x = "Dim1", y = "Dim2", col = "scaled expression") +
      scale_colour_gradient2(low = col_low, mid = col_mid, high = col_high, na.value = "grey", midpoint = 0.5) + 
      ggplot_theme_bw
    }
    return(plot_grid(plot_isoforms@ggplot, plot_expression_umaps, ncol = 1, rel_heights = c(1, 2.5)))
  }
  if (impute) {
    umap_all <- cbind(t(expr_all), umap_all)
    # ggplot2 lazy evaluation with for loops may need to unquote idx with !!
    # operator
    plot_idx <- function(idx) {
      p <- ggplot() + 
      geom_point(data = umap_all, aes(x = x, y = y, col = umap_all[, idx]), alpha = 0.7, size = 0.2) +
      labs(x = "Dim1", y = "Dim2", col = "scaled expression", title = colnames(umap_all)[idx]) +
      scale_colour_gradient2(low = col_low, mid = col_mid, high = col_high, na.value = "grey", 
        midpoint = median(umap_all[, idx], na.rm = TRUE)) +
      ggplot_theme_bw
    }
  } else {
    plot_idx <- function(idx) {
      p <- ggplot() + 
      geom_point(data = umap_80, aes(x = x, y = y), alpha = 0.2, size = 0.2, col = "grey", show.legend = FALSE) +
      geom_point(data = umap_20, aes(x = x, y = y, col = umap_20[, idx]), alpha = 0.7, size = 0.7) +
      labs(x = "Dim1", y = "Dim2", col = "scaled expression", title = colnames(umap_20)[idx]) +
      scale_colour_gradient2(low = col_low, mid = col_mid, high = col_high,
          na.value = "grey", midpoint = median(umap_20[, idx], na.rm = TRUE)) +
      ggplot_theme_bw
      p
    }
  }

  plot_expression_umaps <- lapply(1:n_isoforms, plot_idx)
  plot_expression_umaps <- plot_expression_umaps[match(tr_order, colnames(umap_20)[1:n_isoforms])]

  legend <- ggplot(data.frame(x = 1:2, y = 1:2)) + geom_point(aes(x = x, y = y,
    col = c(-1, 1))) + labs(col = "scaled expression") + scale_colour_gradient2(low = col_low,
    mid = col_mid, high = col_high, breaks = c(-0.8, 0.8), labels = c("Low",
      "High"))
  legend <- get_legend(legend + theme(legend.box.margin = margin(0, 0, 0, 12)))
  combined_isoform_plot <- plot_grid(plotlist = c(list(annotation = plot_isoforms@ggplot),
    plot_expression_umaps))
  combined_isoform_plot <- plot_grid(combined_isoform_plot, legend, rel_widths = c(3,
    0.4))
  return(combined_isoform_plot)
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
#' source(system.file("examples/scmixology1.R", package = "FLAMES"))
#' scm_lib80 <- scuttle::addPerCellQC(scm_lib80)
#' scm_lib20 <- scuttle::addPerCellQC(scm_lib20)
#' qc_80 <- scuttle::quickPerCellQC(colData(scm_lib80))
#' qc_20 <- scuttle::quickPerCellQC(colData(scm_lib20))
#' qc_80$discard[scm_lib80$sum < 10000] <- TRUE
#' qc_20$discard[scm_lib20$sum < 20000] <- TRUE
#' 
#' combined_sce <- combine_sce(
#'     short_read_large = scm_lib80[,!qc_80$discard],
#'     short_read_small = scm_lib20[,!qc_20$discard],
#'     long_read_sce = scm_lib20_transcripts,
#'     remove_duplicates = FALSE)
#'
#' sc_heatmap_expression(gene = "ENSG00000108107", multiAssay = combined_sce)
#' 
#' @export
#' @md
sc_heatmap_expression <- function(gene, multiAssay, impute = FALSE, n_isoforms = 4, transcript_ids, n_pcs = 40,
  isoform_legend_width = 7, col_low = "#313695", col_mid = "#FFFFBF", col_high = "#A50026", color_quantile = 0.95) {

  multiAssay <- sc_reduce_dims(multiAssay, n_pcs)

  sample_lib <- factor(multiAssay[, , "gene_counts"]$Lib_small)
  levels(sample_lib) <- gsub("FALSE", "large", gsub("TRUE", "small", levels(sample_lib)))
  names(sample_lib) <- colnames(multiAssay[, , "gene_counts"])[["gene_counts"]]

  # row_meta columns: 'transcript_id', 'gene_id', 'FSM_match', 'mean',
  # 'detected', 'gene_name' rows: transcript_id
  row_meta <- as.data.frame(rowData(experiments(multiAssay)$transcript_counts))
  if (!missing(transcript_ids) && !is.null(transcript_ids)) {
    if (!all(transcript_ids %in% rownames(row_meta))) {
      message(paste0(c("transcript_id(s) ", paste(transcript_ids[!(transcript_ids %in%
        rownames(transcript_id))], sep = " "), " not found\n")))
    }
    row_meta <- row_meta %>%
      dplyr::filter(transcript_id %in% transcript_ids)
    stopifnot(`no transcript_id match` = dim(row_meta)[1] == 0)
  } else {
    row_meta <- row_meta %>%
      dplyr::filter(gene_id == gene) %>%
      slice_max(n = n_isoforms, order_by = mean)  # keep top n isoforms per gene
    if (dim(row_meta)[1] < 2) {
      stop(paste0("Alternative isoform not found for the given gene: ", gene))
    }
  }

  # apply top n to sce
  tr_sce_multi <- experiments(multiAssay[rownames(row_meta),
        names(sample_lib)[which(sample_lib == "small")],
        "transcript_counts"])$transcript_counts 
  if (n_isoforms > nrow(row_meta)) {
    n_isoforms <- nrow(row_meta)
    message(paste0(c("Only ", n_isoforms, " isoforms found for this gene.\n")))
  }

  # Alignment info
  isoform_sel <- rowRanges(experiments(multiAssay)$transcript_counts[row_meta$FSM_match, ])
  # isoform_sel <- isoform_sel[row_meta$FSM_match] # Order by expression levels
  if (length(isoform_sel) == 2) {
    fill_by_isoform <- c(rep(col_low, length(isoform_sel[[1]])), rep(col_high, length(isoform_sel[[2]])))
    plot_isoforms <- ggbio::autoplot(isoform_sel, label = TRUE, fill = fill_by_isoform) +
      theme_bw() + theme(
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")
      )
  } else {
    plot_isoforms <- ggbio::autoplot(isoform_sel, label = TRUE) +
      theme_bw() + theme(
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")
      )
    # ggbio::geom_alignment does not follow the order in isoform_sel
  }

  # Could try using stat_stepping() to get the order of isoforms
  # https://github.com/lawremi/ggbio/issues/149
  tr_order <- rev(names(isoform_sel)[match(1:n_isoforms, plot_isoforms$layers[[length(plot_isoforms$layers)]]$data$stepping)])
  # isoform_sel <- isoform_sel[tr_order]

  legends_heatmap <- list()
  for (i in names(isoform_sel)) {
    p <- ggplot(isoform_sel[[i]]) +
      geom_alignment(
        label = FALSE, range.geom = "rect",
        gap.geom = "arrow", utr.geom = "rect"
      ) +
      xlim(c(min(min(start(isoform_sel))), max(max(end(isoform_sel))))) +
      theme(
        axis.line = element_blank(), axis.text.x = element_blank(),
        axis.text.y = element_blank(), axis.ticks = element_blank(),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        legend.position = "none", panel.background = element_blank(),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), plot.background = element_blank(),
        plot.margin = margin(20, 0, 20, 0), plot.title = element_text(hjust = 0.5)
      ) +
      labs(title = i)
    legends_heatmap[[i]] <- p@ggplot
  }
  legends_heatmap <- legends_heatmap[tr_order]

  flames_outdir <- experiments(multiAssay)$transcript_counts@metadata$OutputFiles$outdir
  if (!is.null(flames_outdir) && file.exists(file.path(flames_outdir, "cluster_annotation.csv"))) {
    cluster_barcode <- read.csv(file.path(flames_outdir, "cluster_annotation.csv"), stringsAsFactors = FALSE)
    cluster_barcode <- cluster_barcode[, c("barcode_seq", "groups")]
    rownames(cluster_barcode) <- cluster_barcode[, "barcode_seq"]
  } else if (!is.null(multiAssay$cell_type)) {
    cluster_barcode <- data.frame(barcode_seq = names(sample_lib), groups = multiAssay[,names(sample_lib),]$cell_type)
    rownames(cluster_barcode) <- cluster_barcode[, "barcode_seq"]
  }

  group_annotation <- function(cluster_barcode) {
    n <- length(table(cluster_barcode))
    if (n > 11 || n < 2) {
      column_anno <- columnAnnotation(df = cluster_barcode)
    } else if (n == 2) {
      palette <- brewer.pal(n = 3, name = heatmap_annotation_colors)[c(1, 3)]
      names(palette) <- names(table(cluster_barcode))
      column_anno <- columnAnnotation(
        df = cluster_barcode,
        col = list(groups = palette)
      )
    } else {
      palette <- brewer.pal(n, name = heatmap_annotation_colors)
      names(palette) <- names(table(cluster_barcode))
      column_anno <- columnAnnotation(
        df = cluster_barcode,
        col = list(groups = palette)
      )
    }
    return(column_anno)
  }

  # Error in hclust(x) : NA/NaN/Inf in foreign function call (arg 10)
  # removing NaN cells from expr_20 and expr_all
  if (impute) {
    snn <- buildSNNGraph(experiments(multiAssay)$gene_counts, use.dimred = "PCA")
    snn_mat <- as_adjacency_matrix(snn, attr = "weight")
    diag(snn_mat) <- ceiling(max(snn_mat))
    expr_heatmap <- sc_impute_expression(logcounts(tr_sce_multi), snn_mat, sample_lib)
  } else {
    expr_heatmap <- matrix(0, dim(tr_sce_multi)[1], sum(sample_lib == "small"))
    rownames(expr_heatmap) <- rownames(tr_sce_multi)
    colnames(expr_heatmap) <- names(sample_lib[sample_lib == "small"])
    expr_heatmap[, colnames(tr_sce_multi)] <- logcounts(tr_sce_multi)
    expr_heatmap <- scale(expr_heatmap)
    expr_heatmap <- t(apply(expr_heatmap, 1, reduce_quantile))
  }
  expr_heatmap <- expr_heatmap[, !apply(expr_heatmap, 2, function(x) {
    any(is.na(x))
  })]

  cell_order <- stats::hclust(stats::dist(t(expr_heatmap)))$order
  if (exists('cluster_barcode')) {
    cluster_barcode <- cluster_barcode[colnames(expr_heatmap)[cell_order], "groups", drop = FALSE]
  }
  expr_heatmap <- expr_heatmap[tr_order, cell_order]

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
      name = ifelse(impute, "imputed expression", "scaled expression"),
      cluster_rows = FALSE, cluster_columns = FALSE, use_raster = FALSE,
      show_column_names = FALSE, show_row_names = FALSE, 
      # https://www.r-bloggers.com/2017/02/use-switch-instead-of-ifelse-to-return-a-null/
      top_annotation = switch(exists('cluster_barcode'), group_annotation(cluster_barcode)),
      left_annotation = rowAnnotation(isoform = isoform_annotation, annotation_name_rot = 0),
      col = expr_color_mapping(expr_heatmap)
    )
  )
}
