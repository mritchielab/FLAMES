#' FLAMES Annotated Plottings
#'
#' Plot isoform exons alignments for a given gene, along with UMAP showing expression levels.
#'
#'
#' @details
#' This function takes the short-read data (as \code{SingleCellExperiment} objects) from both the smaller
#' and the larger libraries to generate a combined UMAP, the expression levels of isoforms (using long read data)
#' are then overlayed on top of the UMAP. SNN inference based on gene counts were performed to impute isoform expression
#' for cells in the larger library.
#'
#' @param sce_all The \code{SingleCellExperiment} object containing gene counts assay for the combined library.
#' Require \code{sce_all$lib} to specify the library of each cell, either "\code{lib20}" or "\code{lib80}".
#' Alternatively, provide \code{sce_20} \code{sce_80} and the function will create \code{sce_all} using \code{BiocGenerics::cbind}.
#' @param sce_80 The \code{SingleCellExperiment} object containing gene counts assay for the larger (80% ~ 90%) library.
#' @param sce_20 The \code{SingleCellExperiment} object containing gene counts assay for the smaller (10% ~ 20%) library.
#' @param path The path to the folder containing outputs of \code{sc_long_pipeline}.
#' @param gene The gene symbol of interest.
#' @param n_isoforms The number of expressed isoforms to keep.
#' @param n_pcs The number of principal components to generate.
#' @param cluster_annotation Path to the cluster annotation CSV (required for heatmap, if \code{cluster_annotation.csv} is not in \code{path} and \code{sce_all$cell_type} does not exist)
#' @param dup_bc Cell barcodes found both in the larger and smaller library, will be used to filter cells in the long-read
#' data. (Filtering long-read data will be implemented in the main pipeline soon)
#' @param return_sce_all Whether to return the processed \code{SingleCellExperiment} object.
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
#'
#' @importFrom dplyr group_by summarise_at slice_max filter
#' @importFrom tidyr gather pivot_wider
#' @importFrom magrittr "%>%"
#' @importFrom S4Vectors DataFrame
#' @importFrom SingleCellExperiment SingleCellExperiment reducedDimNames logcounts
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
#' @export
#' @md
sc_annotate_umap <- function(gene, path, sce_all = NULL, sce_20 = NULL, sce_80 = NULL, n_isoforms = 4, n_pcs = 40, cluster_annotation, dup_bc = NULL, return_sce_all = TRUE,
                             heatmap_annotation_colors = "BrBG", isoform_legend_width = 7, heatmap_color_quantile = 0.95, col_low = "#313695", col_mid = "#FFFFBF", col_high = "#A50026") {
  ## dup_bc: cell barcodes found in both lib_20 and lib_80
  ## the pipeline should filter out these cells (not implemented yet)
  ## this parameter will need to be removed after implementing filtering in the pipeline

  # Check args
  if (missing(path) || missing(gene)) {
    stop("arguments path and gene are required\n")
  }

  if ((missing(sce_20) || missing(sce_80)) && missing(sce_all)) {
    stop("Please provide sce_all, or both sce_20 and sce_80\n")
  }

  if (!is.null(sce_all) && (is.null(sce_all$lib) || !all(sce_all$lib %in% c("lib20", "lib80")))) {
    stop("Please provide sce_all$lib to specify the library of each cell\n
          Values of sce_all$lib should be either \"lib20\" or \"lib80\".\n")
  }

  outputs <- list()

  # Process short reads SCE
  if (is.null(sce_all)) {
    cat("Creating sce_all from sce_20 and sce_80 ...\n")
    rowData(sce_20) <- rowData(sce_20)[!(names(rowData(sce_20)) %in% c("mean", "detected"))]
    rowData(sce_80) <- rowData(sce_80)[!(names(rowData(sce_80)) %in% c("mean", "detected"))]
    sce_all <- cbind(sce_80, sce_20)
    sce_all$lib <- c(rep("lib80", dim(sce_80)[2]), rep("lib20", dim(sce_20)[2]))
  } else {
    cat("Creating sce_20 and sce_80 from sce_all ... \n")
    sce_20 <- sce_all[, sce_all$lib == "lib20"]
    sce_80 <- sce_all[, sce_all$lib == "lib80"]
  }

  # Process barcodes from short reads
  for (sce in list(sce_20, sce_80, sce_all)) {
    if (is.null((colnames(sce))) && is.null(colData(sce)@listData$Barcode)) {
      stop("Please add barcodes to SCEs, as colnames or colData(sce)@listData$Barcode)")
    } else if (is.null((colnames(sce))) && !is.null(colData(sce)@listData$Barcode)) {
      colnames(sce) <- colData(sce)@listData$Barcode
    }
    if (!all(grepl("-1$", colnames(sce)))) {
      cat("Trimming \"-1\" from colnames\n")
      colnames(sce) <- gsub("-1$", "", colnames(sce))
    }
  }

  # SCE from long reads
  ### Transcrips counts from long read data
  cat("Loading long read data ...\n")
  transcript_count <- read.csv(file.path(path, "transcript_count.csv.gz"), stringsAsFactors = FALSE)
  isoform_FSM_annotation <- read.csv(file.path(path, "isoform_FSM_annotation.csv"), stringsAsFactors = FALSE)
  isoform_gff <- rtracklayer::import.gff3(file.path(path, "isoform_annotated.filtered.gff3"))

  transcript_count <- transcript_count[match(isoform_FSM_annotation$transcript_id, transcript_count$transcript_id), ]
  transcript_count$FSM_match <- isoform_FSM_annotation$FSM_match
  cell_bcs <- colnames(transcript_count)[!(colnames(transcript_count) %in% c("transcript_id", "gene_id", "FSM_match"))]
  tr_anno <- transcript_count[, c("transcript_id", "gene_id", "FSM_match")]

  fsm_csv <- file.path(path, "FSM_count.csv.gz")
  if (file.exists(fsm_csv)) {
    mer_tmp <- read.csv(fsm_csv)
  } else {
    # sum transcript (FSM) counts
    cat("Creating FSM_count.csv.gz ...\n")
    mer_tmp <- transcript_count %>%
      group_by(FSM_match) %>%
      summarise_at(cell_bcs, sum)
    write.csv(mer_tmp, file = gzfile(fsm_csv), row.names = FALSE)
    cat(paste0(c(fsm_csv, "saved.\n")))
  }

  # Create long read SCE
  tr_anno <- tr_anno[match(mer_tmp$FSM_match, tr_anno$FSM_match), ]
  tr_sce <- SingleCellExperiment(assays = list(counts = as.matrix(mer_tmp[, -1])))
  rownames(tr_sce) <- mer_tmp$FSM_match
  rowData(tr_sce) <- DataFrame(tr_anno)

  # QC for long read SCE
  tr_sce <- tr_sce[, !(colnames(tr_sce) %in% dup_bc)]
  tr_sce <- addPerCellQC(tr_sce)
  tr_sce <- addPerFeatureQC(tr_sce)
  keep.hi <- isOutlier(tr_sce$sum, type = "higher", log = TRUE)
  keep.low <- isOutlier(tr_sce$sum, type = "lower", log = TRUE)
  tr_sce <- tr_sce[, (!keep.hi) & (!keep.low)]

  # Handle barcodes
  n_rm <- dim(tr_sce)[2]
  tr_sce <- tr_sce[, colnames(tr_sce) %in% colnames(sce_all)]
  n_rm <- n_rm - dim(tr_sce)[2]
  if (n_rm != 0) {
    message(paste0(c(n_rm, " barcode(s) from long reads are not found in any short read library, they are removed from analysis\n")))
  }

  n_rm <- dim(tr_sce)[2]
  tr_sce <- tr_sce[, colnames(tr_sce) %in% colnames(sce_20)]
  n_rm <- n_rm - dim(tr_sce)[2]
  if (n_rm != 0) {
    message(paste0(c(
      n_rm, " barcode(s) from long reads are assigned to the 80% ~ 90% library (in short read SCEs), they are removed from analysis\n",
      "Please consider revising short read libraries\n"
    )))
  }

  n_rm <- dim(sce_20)[2]
  sce_20 <- sce_20[, colnames(sce_20) %in% colnames(tr_sce)]
  n_rm <- n_rm - dim(sce_20)[2]
  if (n_rm != 0) {
    message(paste0(c(n_rm, " barcode(s) from the 10% ~ 20% library are not found in long reads, removing ... \n")))
    sce_all <- sce_all[, colnames(sce_all) %in% c(colnames(sce_80), colnames(sce_20))]
  }


  ### Process SCEs
  tr_sce <- logNormCounts(tr_sce)
  rowData(tr_sce)$gene_id <- gsub("\\..*", "", rowData(tr_sce)$gene_id)
  # gene_name = mapIds(org.Hs.eg.db,
  #                   keys=rowData(tr_sce)$gene_id,
  #                   column="SYMBOL",
  #                   keytype="ENSEMBL",
  #                   multiVals="first")

  if (!("PCA" %in% reducedDimNames(sce_all)) || dim(sce_all@int_colData$reducedDims$PCA)[2] < n_pcs) {
    cat("Running PCA for sce_all ...\n")
    sce_all <- logNormCounts(sce_all)
    hvgs <- getTopHVGs(modelGeneVar(sce_all), n = 2000)
    sce_all <- fixedPCA(sce_all, rank = n_pcs, subset.row = hvgs)
    snn <- buildSNNGraph(sce_all, use.dimred = "PCA")
  } else if (dim(sce_all@int_colData$reducedDims$PCA)[2] > n_pcs) {
    message(paste0(c(
      "sce_all have ", dim(sce_all@int_colData$reducedDims$PCA)[2], " PCs, using the first ",
      n_pcs, " PCs since n_pcs set to ", n_pcs, "\n"
    )))
    sce_all@int_colData$reducedDims$tmp <- sce_all@int_colData$reducedDims$PCA[, 1:n_pcs]
    snn <- buildSNNGraph(sce_all, use.dimred = "tmp")
    sce_all@int_colData$reducedDims$tmp <- NULL
  } else {
    snn <- buildSNNGraph(sce_all, use.dimred = "PCA")
  }

  # Distance matrix for imputation
  snn_mat <- as_adjacency_matrix(snn, attr = "weight")
  diag(snn_mat) <- ceiling(max(snn_mat))

  if (!("UMAP" %in% reducedDimNames(sce_all))) {
    cat("Running UMAP for sce_all ...\n")
    sce_all <- runUMAP(sce_all, dimred = "PCA")
  } else {
    cat("Skipping runUMAP...\n")
  }

  plot_umap <- ggplot() +
    geom_point(aes(
      x = sce_all@int_colData$reducedDims$UMAP[, 1],
      y = sce_all@int_colData$reducedDims$UMAP[, 2],
      col = factor(sce_all$lib)
    ),
    size = 0.02
    ) +
    labs(x = "Dim1", y = "Dim2", title = "Lib_all UMAP", color = "library")

  #  if ( !("PCA" %in% reducedDimNames(sce_20) )) {
  #    cat("running PCA for sce_20 ...\n")
  #    sce_20 <- logNormCounts(sce_20)
  #    hvgs <- getTopHVGs(modelGeneVar(sce_20), n=2000)
  #    sce_20 <- fixedPCA(sce_20, rank=n_pcs, subset.row=hvgs)
  #  }


  # Select transcript with alternative isoforms
  cat("Finding transcript with alternative isoforms...\n")
  tr_sce_multi <- tr_sce[rowData(tr_sce)$gene_id %in% names(table(rowData(tr_sce)$gene_id)[table(rowData(tr_sce)$gene_id) > 1]), ]

  # row_meta columns: "transcript_id", "gene_id", "FSM_match", "mean", "detected", "gene_name"
  # rows: transcript_id
  row_meta <- as.data.frame(rowData(tr_sce_multi))
  row_meta <- row_meta %>%
    dplyr::filter(gene_id == gene) %>%
    slice_max(n = n_isoforms, order_by = mean) # keep top n isoforms per gene
  if (dim(row_meta)[1] < 2) {
    stop(paste0("Alternative isoform not found for the given gene: ", gene))
  }
  tr_sce_multi <- tr_sce_multi[rowData(tr_sce_multi)$transcript_id %in% row_meta$transcript_id, ] # apply top n to sce
  if (n_isoforms > nrow(row_meta)) {
    n_isoforms <- nrow(row_meta)
    message(paste0(c("Only ", n_isoforms, " isoforms found for this gene.\n")))
  }

  # Alignment info
  isoform_gff$Parent <- as.character(isoform_gff$Parent)
  isoform_gff$transcript_id <- unlist(lapply(strsplit(isoform_gff$Parent, split = ":"), function(x) {
    x[2]
  }))

  isoform_sel <- isoform_gff[isoform_gff$transcript_id %in% rowData(tr_sce_multi)[, "transcript_id"], ]
  isoform_sel <- S4Vectors::split(isoform_sel, isoform_sel$transcript_id) # Split isoforms into a list
  names(isoform_sel) <- row_meta$FSM_match[match(names(isoform_sel), row_meta$transcript_id)] # Set names to FSM
  isoform_sel <- isoform_sel[row_meta$FSM_match] # Order by expression levels
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
  expr <- logcounts(tr_sce_multi)[row_meta$FSM_match, ]
  expr_all <- matrix(0, nrow(expr), ncol(sce_all)) # dim: transcripts x cells, initialised to 0
  rownames(expr_all) <- rownames(expr)
  colnames(expr_all) <- colnames(sce_all)

  # expr <- expr[, colnames(expr) %in% colnames(sce_all)]
  expr_all[, colnames(expr)] <- expr
  expr_all <- expr_all %*% snn_mat
  colnames(expr_all) <- colnames(sce_all)

  expr_all <- t(t(expr_all) / colSums(snn_mat))
  expr_all <- scale(expr_all)

  reduce_quantile <- function(x, q = 0.05) {
    x[x < quantile(x, q, na.rm = TRUE)] <- quantile(x, q, na.rm = TRUE)
    x[x > quantile(x, 1 - q, na.rm = TRUE)] <- quantile(x, 1 - q, na.rm = TRUE)
    x
  }
  expr_all <- t(apply(expr_all, 1, reduce_quantile))

  # Expressions without imputation (Lib_20)
  expr_20 <- matrix(0, nrow(expr), ncol(sce_20))
  rownames(expr_20) <- rownames(expr)
  colnames(expr_20) <- colnames(sce_20)
  expr_20[, colnames(expr)] <- expr
  expr_20 <- scale(expr_20)
  expr_20 <- t(apply(expr_20, 1, reduce_quantile))

  # scale() will produce NaN values for cell with 0 counts for all isoforms
  # using na.value = "grey" for UMAPs

  ### Plot with UMAP from lib_all
  cat("Plotting expression UMAPs ...\n")
  umap_all <- as.data.frame(sce_all@int_colData$reducedDims$UMAP)
  colnames(umap_all) <- c("x", "y")

  umap_80 <- umap_all[sce_all$lib == "lib80", ]
  colnames(umap_80) <- c("x", "y")

  umap_20 <- umap_all[sce_all$lib == "lib20", ]
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
  if (file.exists(file.path(path, "cluster_annotation.csv")) || (!missing(cluster_annotation) && file.exists(cluster_annotation)) || !is.null(sce_all$cell_type)) {
    cat("Plotting heatmaps ...\n")
    if (file.exists(file.path(path, "cluster_annotation.csv"))) {
      cluster_barcode <- read.csv(file.path(path, "cluster_annotation.csv"), stringsAsFactors = FALSE)
      cluster_barcode <- cluster_barcode[, c("barcode_seq", "groups")]
    } else if (!missing(cluster_annotation) && file.exists(cluster_annotation)) {
      cluster_barcode <- read.csv(cluster_annotation, stringsAsFactors = FALSE)
      cluster_barcode <- cluster_barcode[, c("barcode_seq", "groups")]
    } else {
      cluster_barcode <- data.frame(barcode_seq = colnames(sce_all), groups = sce_all$cell_type)
    }
    rownames(cluster_barcode) <- cluster_barcode[, "barcode_seq"]

    outputs[["umap_clusters"]] <- ggplot(umap_all) +
      geom_point(aes(
        x = x, y = y,
        col = factor(cluster_barcode[colnames(sce_all), "groups"])
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
        grid.arrange(grobs = legends_heatmap, ncol = 1, vp = viewport(), newpage = FALSE)
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


  if (return_sce_all) {
    outputs[["sce_all"]] <- sce_all
  }

  return(outputs)
}