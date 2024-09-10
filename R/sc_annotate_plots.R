#' Combine SCE
#'
#' Combine FLT-seq SingleCellExperiment objects
#'
#' @details For protcols like FLT-seq that generate two libraries, one with both short and long reads,
#' and one with only short reads, this function combines the two libraries into a single
#' \code{SingleCellExperiment} object. For the library with both long and short reads, the long-read
#' transcript counts should be stored in the 'transcript' altExp slot of the \code{SingleCellExperiment}
#' object. This function will combine the short-read gene counts of both libraries, and for the
#' transcripts counts, it will leave \code{NA} values for the cells from the short-read only library.
#' The \code{sc_impute_transcript} function can then be used to impute the \code{NA} values.
#'
#' @param sce_with_lr A \code{SingleCellExperiment} object with both long and short reads. The long-read
#' transcript counts should be stored in the 'transcript' altExp slot.
#' @param sce_without_lr A \code{SingleCellExperiment} object with only short reads.
#'
#' @return A \code{SingleCellExperiment} object with combined gene counts and a "transcript" altExp slot.
#'
#' @examples
#' with_lr <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = matrix(rpois(100, 5), ncol = 10)))
#' without_lr <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = matrix(rpois(200, 5), ncol = 20)))
#' long_read <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = matrix(rpois(50, 5), ncol = 10)))
#' SingleCellExperiment::altExp(with_lr, "transcript") <- long_read
#' SummarizedExperiment::colData(with_lr)$Barcode <- paste0(1:10, "-1")
#' SummarizedExperiment::colData(without_lr)$Barcode <- paste0(8:27, "-1")
#' rownames(with_lr) <- as.character(101:110)
#' rownames(without_lr) <- as.character(103:112)
#' rownames(long_read) <- as.character(1001:1005)
#' combined_sce <- FLAMES::combine_sce(sce_with_lr = with_lr, sce_without_lr = without_lr)
#' combined_sce
#'
#' @importFrom BiocGenerics cbind colnames rownames
#' @importFrom SingleCellExperiment SingleCellExperiment altExp altExp<- counts counts<-
#' @importFrom SummarizedExperiment rowData colData rowData<- colData<- rowRanges rowRanges<-
#' @importFrom DropletUtils read10xCounts
#' @export
#' @md
combine_sce <- function(sce_with_lr, sce_without_lr) {
  if (is.character(sce_without_lr)) {
    sce_without_lr <- read10xCounts(sce_without_lr)
  }

  if (!is.null(colnames(sce_with_lr)) & !is.null(colnames(sce_without_lr))) {
    duplicated_barcodes <- intersect(colnames(sce_with_lr), colnames(sce_without_lr))
    sce_without_lr <- sce_without_lr[, !colnames(sce_without_lr) %in% duplicated_barcodes]
    newColNames <- c(colnames(sce_with_lr), colnames(sce_without_lr))
  } else if (!is.null(colData(sce_with_lr)$Barcode) & !is.null(colData(sce_without_lr)$Barcode)) {
    duplicated_barcodes <- intersect(colData(sce_with_lr)$Barcode, colData(sce_without_lr)$Barcode)
    sce_without_lr <- sce_without_lr[, !colData(sce_without_lr)$Barcode %in% duplicated_barcodes]
    newColNames <- c(colData(sce_with_lr)$Barcode, colData(sce_without_lr)$Barcode) 
    message("Missing colnames, using Barcode from colData as colnames")
  } else {
    stop("Barcode information not found, please either provide it as colData(sce)$Barcode or colnames(sce)")
  }

  gene_intersect <- intersect(rownames(sce_with_lr), rownames(sce_without_lr))
  sce_with_lr <- sce_with_lr[gene_intersect, ]
  sce_without_lr <- sce_without_lr[gene_intersect, ]

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

  colnames(combined_sce) <- newColNames
  tryCatch({
    colDataCols <- intersect(colnames(colData(sce_with_lr)), colnames(colData(sce_without_lr)))
    colData(combined_sce) <- rbind(
      colData(sce_with_lr)[,colDataCols, drop=FALSE], 
      colData(sce_without_lr)[,colDataCols, drop=FALSE]
    )
  }, error = function(e) {
    message("Failed to combine colData, error: ", e)
    message("Returning combined SCE without colData")
  })
  combined_sce
}

#' Impute missing transcript counts
#'
#' Impute missing transcript counts using a shared nearest neighbor graph
#'
#' @details For cells with \code{NA} values in the "transcript" altExp slot, this function imputes the
#' missing values from cells with non-missing values. A shared nearest neighbor graph is built using 
#' reduced dimensions from the \code{SingleCellExperiment} object, and the imputation is done where
#' the imputed value for a cell is the weighted sum of the transcript counts of its neighbors.
#' Imputed values are stored in the "logcounts" assay of the "transcript" altExp slot.
#' The "counts" assay is used to obtain logcounts but left unchanged.
#'
#' @param combined_sce A \code{SingleCellExperiment} object with gene counts and a "transcript" altExp slot.
#' @param dimred The name of the reduced dimension to use for building the shared nearest neighbor graph.
#' @param ... Additional arguments to pass to \code{scran::buildSNNGraph}. E.g. \code{k = 30}.
#'
#' @return A \code{SingleCellExperiment} object with imputed logcounts assay in the "transcript" altExp slot.
#'
#' @examples
#' sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = matrix(rpois(50, 5), ncol = 10)))
#' long_read <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = matrix(rpois(40, 5), ncol = 10)))
#' SingleCellExperiment::altExp(sce, "transcript") <- long_read
#' SingleCellExperiment::counts(SingleCellExperiment::altExp(sce))[,1:2] <- NA
#' SingleCellExperiment::counts(SingleCellExperiment::altExp(sce))
#' imputed_sce <- sc_impute_transcript(sce, k = 4)
#' SingleCellExperiment::logcounts(SingleCellExperiment::altExp(imputed_sce))
#'
#' @importFrom SingleCellExperiment altExp altExp<- counts counts<- altExps reducedDimNames
#' @importFrom SummarizedExperiment assays assay<-
#' @importFrom scran buildSNNGraph
#' @importFrom scater runPCA runUMAP
#' @importFrom scuttle logNormCounts
#' @importFrom igraph as_adjacency_matrix
#' @export
#' @md
sc_impute_transcript <- function(combined_sce, dimred = "PCA", ...) {
  if (!"transcript" %in% names(SingleCellExperiment::altExps(combined_sce))) {
    stop("transcript counts not found")
  }
  if (!"logcounts" %in% names(SummarizedExperiment::assays(combined_sce))) {
    combined_sce <- scuttle::logNormCounts(combined_sce)
  }
  if (!"PCA" %in% SingleCellExperiment::reducedDimNames(combined_sce)) {
    combined_sce <- scater::runPCA(combined_sce)
  }
  if (!"UMAP" %in% SingleCellExperiment::reducedDimNames(combined_sce) & dimred == "UMAP") {
    combined_sce <- scater::runUMAP(combined_sce)
  }
  snn <- scran::buildSNNGraph(combined_sce,use.dimred = dimred, ...)
  snn_mat <- igraph::as_adjacency_matrix(snn, attr = "weight")

  imputed_count <- impute(counts(altExp(combined_sce, "transcript")), snn_mat)
  assay(altExp(combined_sce, "transcript"), "logcounts") <- imputed_count
  combined_sce
}

#' @importFrom scuttle normalizeCounts
#' @importFrom BiocGenerics t
impute <- function(mtx, distance_matrix) {
  cat("Imputing transcript counts ...\n")
  expr_impute <- matrix(0, nrow(mtx), ncol(mtx))
  rownames(expr_impute) <- rownames(mtx)
  colnames(expr_impute) <- colnames(mtx)

  cols_with_counts <- apply(mtx, 2, function(x) {
    !any(is.na(x))
  })
  
  real_counts <- scuttle::normalizeCounts(mtx[, cols_with_counts])
  expr_impute[, cols_with_counts] <- real_counts

  expr_impute <- expr_impute %*% distance_matrix
  colnames(expr_impute) <- colnames(mtx)
  expr_impute <- t(t(expr_impute) / (distance_matrix%*%cols_with_counts))

  expr_impute[, cols_with_counts] <- real_counts
  if (any(is.infinite(expr_impute))) {
    # in case of division by 0
    warning("Failed to impute some cells, leaving them as NA")
    expr_impute[is.infinite(expr_impute)] <- NA
  }

  expr_impute
}


#' @importFrom stats quantile
reduce_quantile <- function(x, q = 0.05) {
  x[x < quantile(x, q, na.rm = TRUE)] <- quantile(x, q, na.rm = TRUE)
  x[x > quantile(x, 1 - q, na.rm = TRUE)] <- quantile(x, 1 - q, na.rm = TRUE)
  x
}

# @importFrom SingleCellExperiment counts rowData
get_top_transcript_ids <- function(sce, gene_id, transcript_ids, n) {
  if (missing(transcript_ids)) {
    sce <- sce[rowData(sce)$gene_id == gene_id, ]
    sce <- sce[order(rowSums(counts(sce)), decreasing = TRUE), ]
    if (nrow(sce) > n) {
      sce <- sce[1:n, ]
    }
  } else {
    return(transcript_ids)
  }
  return(rowData(sce)$transcript_id)
}

#' Plot isoforms
#'
#' Plot isoforms, either from a gene or a list of transcript ids.
#'
#' @details
#' This function takes a \code{SingleCellExperiment} object and plots the top isoforms of a gene,
#' or a list of specified transcript ids. Either as a list of plots or together in a grid.
#' This function wraps the \code{ggbio::geom_alignment} function to plot the isoforms, and orders
#' the isoforms by expression levels (when specifying a gene) or by the order of the transcript_ids.
#'
#' @param sce The \code{SingleCellExperiment} object containing transcript counts, 
#' \code{rowRanges} and \code{rowData} with \code{gene_id} and \code{transcript_id} columns.
#' @param gene_id The gene symbol of interest, ignored if \code{transcript_ids} is provided.
#' @param transcript_ids The transcript ids to plot.
#' @param n The number of top isoforms to plot from the gene. Ignored if \code{transcript_ids} is provided.
#' @param format The format of the output, either "plot_grid" or "list".
#'
#' @return When \code{format = "list"}, a list of \code{ggplot} objects is returned. 
#' Otherwise, a grid of the plots is returned.
#'
#' @examples
#' plot_isoforms(scmixology_lib10_transcripts, gene_id = "ENSG00000108107")
#' 
#' @importFrom SummarizedExperiment rowRanges rowData
#' @importFrom ggbio geom_alignment
#' @importFrom ggplot2 ggplot aes theme_void theme element_line element_text xlim
#' @importFrom cowplot plot_grid
#' @importFrom GenomeInfoDb seqnames
#' @export
plot_isoforms <- function(sce, gene_id, transcript_ids, n = 4, format = "plot_grid") {
  transcript_ids <- get_top_transcript_ids(sce, gene_id, transcript_ids, n)
  sce <- sce[match(transcript_ids, rowData(sce)$transcript_id), ]
  transcripts <- rowRanges(sce)
  x_range <- c(
    min(min(start(transcripts))),
    max(max(end(transcripts)))
  )

  chr <- lapply(seqnames(transcripts), function(x) {
    as.character(unique(x))
  }) |>
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
  if (format == "list") {
    return(plot_list)
  }
  return(cowplot::plot_grid(plotlist = plot_list, ncol = 1))
}

#' FLAMES heetmap plots
#'
#' Plot expression heatmap of top n isoforms of a gene
#'
#' @details
#' Takes \code{SingleCellExperiment} object and plots an expression heatmap with the 
#' isoform visualizations along genomic coordinates.
#'
#' @param sce The \code{SingleCellExperiment} object containing transcript counts, 
#' \code{rowRanges} and \code{rowData} with \code{gene_id} and \code{transcript_id} columns.
#' @param gene_id The gene symbol of interest, ignored if \code{transcript_ids} is provided.
#' @param transcript_ids The transcript ids to plot.
#' @param n The number of top isoforms to plot from the gene. Ignored if \code{transcript_ids} is provided.
#' @param isoform_legend_width The width of isoform legends in heatmaps, in \code{cm}.
#' @param col_low Color for cells with low expression levels in UMAPs.
#' @param col_mid Color for cells with intermediate expression levels in UMAPs.
#' @param col_high Color for cells with high expression levels in UMAPs.
#' @param color_quantile The lower and upper expression quantile to be displayed bewteen \code{col_low} and \code{col_high}, e.g. with \code{color_quantile = 0.95}, cells with expressions higher than 95% of other cells will all be shown in \code{col_high}, and cells with expression lower than 95% of other cells will all be shown in \code{col_low}.
#'
#' @return a \code{ComplexHeatmap}
#'
#' @examples
#' scmixology_lib10_transcripts |>
#'   scuttle::logNormCounts() |>
#'   plot_isoform_heatmap(gene = "ENSG00000108107")
#'
#' @importFrom SingleCellExperiment rowData logcounts colLabels
#' @importFrom ComplexHeatmap AnnotationFunction Heatmap rowAnnotation
#' @importFrom grid unit viewport
#' @importFrom gridExtra grid.arrange
#' @importFrom RColorBrewer brewer.pal
#' @importFrom circlize colorRamp2
#' @importFrom stats quantile hclust dist
#'
#' @export
#' @md
plot_isoform_heatmap <- function(
    sce, gene_id, transcript_ids, n = 4,
    isoform_legend_width = 7, col_low = "#313695", col_mid = "#FFFFBF", col_high = "#A50026", color_quantile = 0.95) {
  transcript_ids <- get_top_transcript_ids(sce, gene_id, transcript_ids, n)
  sce <- sce[match(transcript_ids, rowData(sce)$transcript_id), ]
  legends_heatmap <- plot_isoforms(sce, gene_id, transcript_ids, n, format = "list")

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

  expr_color_mapping <- function(expr_matrix) {
    if (color_quantile > 1 || color_quantile < 0) {
      color_quantile <- 0.95
    }
    breaks <- rep(0, 3)
    breaks[1] <- quantile(expr_matrix, 1- color_quantile, na.rm = TRUE)
    breaks[3] <- quantile(expr_matrix, color_quantile, na.rm = TRUE)
    breaks[2] <- mean(breaks[c(1, 3)])
    return(colorRamp2(breaks, c(col_low, col_mid, col_high)))
  }

  isoform_annotation <- AnnotationFunction(
    fun = function(index) {
      gridExtra::grid.arrange(
        grobs = legends_heatmap, ncol = 1, vp = viewport(),
        newpage = FALSE
      )
    }, which = "row", width = unit(isoform_legend_width, "cm"), n = length(legends_heatmap),
    subsettable = FALSE
  )

  return(
    Heatmap(logcounts(sce),
      name = "log expression",
      cluster_rows = FALSE, cluster_columns = FALSE, use_raster = FALSE,
      show_column_names = FALSE, show_row_names = FALSE,
      # https://www.r-bloggers.com/2017/02/use-switch-instead-of-ifelse-to-return-a-null/
      top_annotation = switch(!is.null(colLabels(sce)),
        group_annotation(colLabels(sce))
      ),
      left_annotation = rowAnnotation(isoform = isoform_annotation, annotation_name_rot = 0),
      col = expr_color_mapping(logcounts(sce))
    )
  )
}

#' FLAMES isoform reduced dimensions plots
#'
#' Plot expression of top n isoforms of a gene in reduced dimensions
#'
#' @details
#' Takes \code{SingleCellExperiment} object and plots an expression on reduced dimensions 
#' with the isoform visualizations along genomic coordinates.
#'
#' @param sce The \code{SingleCellExperiment} object containing transcript counts, 
#' \code{rowRanges} and \code{rowData} with \code{gene_id} and \code{transcript_id} columns.
#' @param gene_id The gene symbol of interest, ignored if \code{transcript_ids} is provided.
#' @param transcript_ids The transcript ids to plot.
#' @param n The number of top isoforms to plot from the gene. Ignored if \code{transcript_ids} is provided.
#' @param reduced_dim_name The name of the reduced dimension to use for plotting cells.
#' @param use_gene_dimred Whether to use gene-level reduced dimensions for plotting. Set to \code{TRUE} 
#' if the \code{SingleCellExperiment} has gene counts in main assay and transcript counts in \code{altExp}.
#' @param format The format of the output, either "plot_grid" or "list".
#' @param expr_func The function to extract expression values from the \code{SingleCellExperiment} object.
#' Default is \code{logcounts}. Alternatively, \code{counts} can be used for raw counts.
#' @param col_low Color for cells with low expression levels in UMAPs.
#' @param col_mid Color for cells with intermediate expression levels in UMAPs.
#' @param col_high Color for cells with high expression levels in UMAPs.
#' @param ... Additional arguments to pass to \code{plot_grid}.
#'
#' @return a \code{ggplot} object of the UMAP(s)
#'
#' @examples
#' scmixology_lib10 <- 
#'   scmixology_lib10[, colSums(SingleCellExperiment::counts(scmixology_lib10)) > 0]
#' sce_lr <- scmixology_lib10[, colnames(scmixology_lib10) %in% colnames(scmixology_lib10_transcripts)]
#' SingleCellExperiment::altExp(sce_lr, "transcript") <-
#'   scmixology_lib10_transcripts[, colnames(sce_lr)]
#' combined_sce <- combine_sce(sce_lr, scmixology_lib90)
#' combined_sce <- combined_sce |>
#'   scuttle::logNormCounts() |>
#'   scater::runPCA() |>
#'   scater::runUMAP()
#' combined_imputed_sce <- sc_impute_transcript(combined_sce)
#' plot_isoform_reduced_dim(combined_sce, 'ENSG00000108107')
#' plot_isoform_reduced_dim(combined_imputed_sce, 'ENSG00000108107')
#'
#' @importFrom SingleCellExperiment logcounts altExpNames altExp reducedDim counts
#' @importFrom SummarizedExperiment rowData assayNames
#' @importFrom BiocGenerics colnames rownames
#' @importFrom scuttle normalizeCounts
#' @importFrom cowplot plot_grid
#' @importFrom ggplot2 aes element_line element_text
#' @importFrom ggplot2 ggplot geom_point labs scale_colour_gradient2 theme_bw
#' @export
#' @md
plot_isoform_reduced_dim <- function(
    sce, gene_id, transcript_ids, n = 4, reduced_dim_name = "UMAP",
    use_gene_dimred = FALSE, expr_func = function(x) {
      SingleCellExperiment::logcounts(x)
    },
    col_low = "#313695", col_mid = "#FFFFBF", col_high = "#A50026", format = "plot_grid", ...) {

  if (!"transcript_id" %in% colnames(rowData(sce)) & "transcript" %in% altExpNames(sce)) {
    use_gene_dimred <- TRUE
  }

  df <- data.frame(
    x = reducedDim(sce, reduced_dim_name)[, 1],
    y = reducedDim(sce, reduced_dim_name)[, 2]
  )

  if (use_gene_dimred) {
    sce <- altExp(sce, "transcript")
    transcript_ids <- get_top_transcript_ids(sce, gene_id, transcript_ids, n)
    sce <- sce[match(transcript_ids, rowData(sce)$transcript_id), ]
    if (!"logcounts" %in% assayNames(sce)) {
      # normalize counts if not already done
      # handle NAs
      mtx <- matrix(NA, nrow(sce), ncol(sce))
      cols_with_counts <- apply(counts(sce), 2, function(x) {
        !any(is.na(x))
      })
      cols_with_counts <- cols_with_counts & (colSums(counts(sce) > 0))
      mtx[, cols_with_counts] <- scuttle::normalizeCounts(counts(sce)[, cols_with_counts])
      colnames(mtx) <- colnames(sce)
      rownames(mtx) <- rownames(sce)
      SingleCellExperiment::logcounts(sce) <- mtx
    }
  } else {
    transcript_ids <- get_top_transcript_ids(sce, gene_id, transcript_ids, n)
    sce <- sce[match(transcript_ids, rowData(sce)$transcript_id), ]
  }
  isoform_plot <- plot_isoforms(sce, gene_id, transcript_ids, n)

  umaps <- list()
  for (i in seq_along(transcript_ids)) {
    df$expr <- expr_func(sce)[i, ]
    p <- ggplot(df) +
      geom_point(aes(x = x, y = y, col = expr), alpha = 0.7, size = 0.2) +
      labs(x = "Dim1", y = "Dim2", col = "expression") +
      scale_colour_gradient2(
        low = col_low, mid = col_mid, high = col_high,
        na.value = "grey",  limits = quantile(df$expr, na.rm = TRUE, c(0.05, 0.95)),
        midpoint = mean(quantile(df$expr, na.rm = TRUE, c(0.05, 0.95))),
      ) +
      theme_bw()
    umaps <- append(umaps, list(p))
  }

  plot_list <- c(list(isoform_plot), umaps)
  if (format == "list") {
    return(umaps)
  }
  return(cowplot::plot_grid(plotlist = plot_list, ...))
}
