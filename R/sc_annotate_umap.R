# TODO:
#       Better visualisation
#       add documentation
#
#' FLAMES Annotated Plottings
#' 
#' Plot isoform exons alignments for a given gene, along with UMAP showing expression levels.
#' 
#' 
#' @details
#' This function takes the short-read data (as \code{SingleCellExperiment} objects) from both the smaller
#' and the larger libraries to generate a combined UMAP, the expression levels of isoforms (using long read data)
#' are then overlayed on top of the UMAP. SNN inference based on gene counts were performed to impute isoform expression
#' for cells in the smaller library.
#' 
#' @param sce_80 The \code{SingleCellExperiment} object containing gene counts assay for the larger library.
#' @param sce_20 The \code{SingleCellExperiment} object containing gene counts assay for the smaller library.
#' @param sce_all The \code{SingleCellExperiment} object containing gene counts assay for the combined library.
#' @param path The path to the folder containing outputs of \code{sc_long_pipeline}.
#' @param gene The gene symbol of interest.
#' @param n_isoforms The number of expressed isoforms to keep.
#' @param n_pcs The number of principal components to generate.
#' @param cluster_annotation Path to the cluster annotation CSV (required for heatmap, if \code{cluster_annotation.csv} is not in \code{path} and \code{sce_20$cell_type} does not exist)
#' @param dup_bc Cell barcodes found both in the larger and smaller library, will be used to filter cells in the long-read
#' data. (Filtering long-read data will be implemented in the main pipeline soon)
#' @param return_sce_all Whether to return the processed \code{SingleCellExperiment} object.
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
#' @importFrom ggplot2 ggplot geom_point aes labs element_blank element_line theme_bw theme scale_colour_gradient2 margin
#' @importFrom ggbio autoplot
#' @importFrom GenomicRanges GRanges GRangesList
#' @importFrom rtracklayer import.gff3
#' @importFrom BiocGenerics cbind colnames rownames
#' @importFrom pheatmap pheatmap
#' @importFrom igraph as_adjacency_matrix
#' @importFrom rtracklayer import
#' @importFrom Matrix t colSums
#' @importFrom cowplot plot_grid get_legend
#' @importFrom ggplotify as.ggplot
#' @export
sc_annotate_umap <- function(sce_20, sce_80, path, gene, n_isoforms = 4, n_pcs = 20, cluster_annotation, dup_bc=NULL, sce_all = NULL, return_sce_all = T){
  ## dup_bc: cell barcodes found in both lib_20 and lib_80
  ## the pipeline should filter out these cells (not implemented yet)
  ## this parameter will need to be removed after implementing filtering in the pipeline

  if (any(missing(sce_20), missing(path), missing(gene))){
    stop("arguments sce_20, path and gene are required\n")
  }

  if (missing(sce_80) && is.null(sce_all)) {
    stop("please provide at least one of:\n - sce_80, or\n - sce_all\n")
  }

  if (is.null(sce_all)) {
    cat("creating sce_all from sce_20 and sce_80 ...\n")
    rowData(sce_20) <- rowData(sce_20)[!(names(rowData(sce_20)) %in% c("mean","detected"))]
    rowData(sce_80) <- rowData(sce_80)[!(names(rowData(sce_80)) %in% c("mean","detected"))]
    sce_all<-cbind(sce_80, sce_20)
    sce_all$lib <- c(rep('lib80', dim(sce_80)[2]), rep('lib20', dim(sce_20)[2]))
  } else if (is.null(sce_all$lib)) {
    sce_all$lib <- rep("lib80", dim(sce_all)[2])
    sce_all$lib[match(rownames(sce_20), rownames(sce_all))] <- "lib20"
  }

  if ( !("PCA" %in% reducedDimNames(sce_all)) && !("UMAP" %in% reducedDimNames(sce_all))) {
    cat("Running PCA for sce_all ...\n")
    sce_all <- logNormCounts(sce_all)
    hvgs <- getTopHVGs(modelGeneVar(sce_all), n=2000)
    sce_all <- fixedPCA(sce_all, rank=n_pcs, subset.row=hvgs)
  }

  if ( !("UMAP" %in% reducedDimNames(sce_all)) ) {
    cat("Running UMAP for sce_all ...\n")
    sce_all <- runUMAP(sce_all)
  }

  plot_umap <- ggplot()+
    geom_point(aes(x=sce_all@int_colData$reducedDims$UMAP[,1], 
                   y=sce_all@int_colData$reducedDims$UMAP[,2], 
                   col=factor(sce_all$lib)),
               size=0.02)+
    labs(x="Dim1", y="Dim2", title = 'Lib_all UMAP', color = "library")
 
  if ( !("PCA" %in% reducedDimNames(sce_20) )) {
    cat("running PCA for sce_20 ...\n")
    sce_20 <- logNormCounts(sce_20)
    hvgs <- getTopHVGs(modelGeneVar(sce_20), n=2000)
    sce_20 <- fixedPCA(sce_20, rank=n_pcs, subset.row=hvgs)
  }
  snn <- buildSNNGraph(sce_20, use.dimred = "PCA")
  snn_mat <- as_adjacency_matrix(snn, attr="weight")
  diag(snn_mat) <- ceiling(max(snn_mat))
 
  ### Transcrips counts from long read data
  cat("Loading long read data ...\n")
  transcript_count <- read.csv(file.path(path,"transcript_count.csv.gz"), stringsAsFactors=FALSE)
  isoform_FSM_annotation <- read.csv(file.path(path,"isoform_FSM_annotation.csv") , stringsAsFactors=FALSE)
  isoform_gff <- rtracklayer::import.gff3(file.path(path,"isoform_annotated.filtered.gff3"))
 
  transcript_count <- transcript_count[match(isoform_FSM_annotation$transcript_id,transcript_count$transcript_id),]
  transcript_count$FSM_match <- isoform_FSM_annotation$FSM_match
  cell_bcs <- colnames(transcript_count)[!(colnames(transcript_count) %in% c("transcript_id","gene_id","FSM_match"))]
  tr_anno <- transcript_count[,c("transcript_id","gene_id","FSM_match")]

  fsm_csv <- file.path(path, "FSM_count.csv.gz")
  if (file.exists(fsm_csv)) {
    mer_tmp <- read.csv(fsm_csv)
  } else {
    #sum transcript (FSM) counts
    cat("Creating FSM_count.csv.gz ...\n")
    mer_tmp <- transcript_count %>%
      group_by(FSM_match) %>%
      summarise_at(cell_bcs,sum)
    write.csv(mer_tmp, file = gzfile(fsm_csv), row.names=FALSE)
    cat(paste0(c(fsm_csv, "saved.\n")))
  }

 
  tr_anno <- tr_anno[match(mer_tmp$FSM_match,tr_anno$FSM_match),]
  tr_sce <- SingleCellExperiment(assays=list(counts=as.matrix(mer_tmp[,-1]) ))
  rownames(tr_sce) <- mer_tmp$FSM_match
  rowData(tr_sce) <- DataFrame(tr_anno)
  tr_sce <- tr_sce[,!(colnames(tr_sce) %in% dup_bc)]
  tr_sce <- addPerCellQC(tr_sce)
  tr_sce <- addPerFeatureQC(tr_sce)
  keep.hi <- isOutlier(tr_sce$sum, type="higher", log=TRUE)
  keep.low <- isOutlier(tr_sce$sum, type="lower", log=TRUE)
  tr_sce <- tr_sce[,(!keep.hi) & (!keep.low)]
  tr_sce <- logNormCounts(tr_sce)
 
  rowData(tr_sce)$gene_id <- gsub("\\..*","", rowData(tr_sce)$gene_id)
  #gene_name = mapIds(org.Hs.eg.db,
  #                   keys=rowData(tr_sce)$gene_id,
  #                   column="SYMBOL",
  #                   keytype="ENSEMBL",
  #                   multiVals="first")
 
  # Select transcript with alternative isoforms
  cat("Finding transcript with alternative isoforms...\n")
  tr_sce_multi <- tr_sce[rowData(tr_sce)$gene_id %in% names( table(rowData(tr_sce)$gene_id)[ table(rowData(tr_sce)$gene_id)>1]), ]
 
  # row_meta columns: "transcript_id", "gene_id", "FSM_match", "mean", "detected", "gene_name"
  # rows: transcript_id
  row_meta <- as.data.frame(rowData(tr_sce_multi))
  row_meta <- row_meta %>% dplyr::filter(gene_id == gene) %>% slice_max(n = n_isoforms, order_by = mean) #keep top n isoforms per gene
  if (dim(row_meta)[1] < 2) {
    stop(paste0("Alternative isoform not found for the given gene: ", gene))
  } 
  tr_sce_multi <- tr_sce_multi[rowData(tr_sce_multi)$transcript_id %in% row_meta$transcript_id,] # apply top n to sce
  fsm_ids <- row_meta$FSM_match 
  if (n_isoforms > length(fsm_ids)) {
    n_isoforms <- length(fsm_ids)
  }
  
  isoform_gff$Parent <- as.character(isoform_gff$Parent)
  isoform_gff$transcript_id <- unlist(lapply(strsplit(isoform_gff$Parent, split = ":"),function(x){x[2]}))
  tr_ids <- rowData(tr_sce_multi)[,"transcript_id"] 

  isoform_sel <- isoform_gff[isoform_gff$transcript_id %in% tr_ids,] 
  isoform_sel <- S4Vectors::split(isoform_sel, isoform_sel$transcript_id) # Split isoforms into a list
  names(isoform_sel) <- fsm_ids[match(names(isoform_sel),tr_ids)] #Set names to FSM
  isoform_sel <- isoform_sel[fsm_ids] # Order by expression levels
  if (length(isoform_sel) == 2) {
    fill_by_isoform <- c( rep("#A50026", length(isoform_sel[[1]])), rep("#313695", length(isoform_sel[[2]])) )
    plot_isoforms <- ggbio::autoplot(isoform_sel, label = TRUE, fill = fill_by_isoform) + 
      theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                       panel.background = element_blank(), axis.line = element_line(colour = "black"))
    plot_isoforms_plain <- ggbio::autoplot(isoform_sel, label = TRUE) + 
      theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                       panel.background = element_blank(), axis.line = element_line(colour = "black"))
  } else {
    plot_isoforms <- ggbio::autoplot(isoform_sel, label = TRUE) + 
      theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                       panel.background = element_blank(), axis.line = element_line(colour = "black"))
    # ggbio::geom_alignment does not follow the order in isoform_sel ...
  }

  # ggbio::geom_alignment orders isoforms by sum of exon lengths
  tr_len <- rep(0,length(isoform_sel))
  for (i in 1:length(isoform_sel)) {
    tr_len[i] <- sum(isoform_sel[[i]]@ranges@width)
  }
  tr_order <- names(isoform_sel)[match(sort(tr_len, decreasing = F), tr_len)]

 
  ### impute transcript counts for cells in lib_20
  cat("Imputing transcript counts for cells in lib_20 ...\n")
  expr <- logcounts(tr_sce_multi)[fsm_ids,]
  expr_20 <- matrix(0,nrow(expr),ncol(sce_20))  #dim: transcripts x cells, initialised to 0
  rownames(expr_20) <- rownames(expr)
  colnames(expr_20) <- colnames(sce_20)
 
  #unsampled <- colnames(expr)[!(colnames(expr) %in% colnames(sce_20))]
  expr <- expr[,colnames(expr) %in% colnames(sce_20)]
  expr_20[,colnames(expr)] <- expr
  expr_20 <- expr_20 %*% snn_mat
  colnames(expr_20) <- colnames(sce_20)
  
  expr_20 <- t(t(expr_20)/colSums(snn_mat))
  expr_20 <-  scale(expr_20)

  reduce_quantile <- function(x, q=0.05) {
    x[x<quantile(x, q, na.rm=T)] <- quantile(x, q, na.rm=T)
    x[x>quantile(x, 1-q, na.rm=T)] <- quantile(x, 1-q, na.rm=T)
    x
  }
  expr_20 <- t(apply(expr_20,1, reduce_quantile))
 
 
  ### Plot with UMAP from lib_all
  cat("Plotting expression UMAPs")
  umap_80 <- as.data.frame(sce_all@int_colData$reducedDims$UMAP)
  umap_80 <- umap_80[sce_all$lib == "lib80",]
  colnames(umap_80) <- c("x", "y")

  umap_20 <- as.data.frame(sce_all@int_colData$reducedDims$UMAP)
  umap_20 <- umap_20[sce_all$lib == "lib20",]
  colnames(umap_20) <- c("x", "y") 

  if (n_isoforms == 2) {
    umap_20 <- cbind(expr=expr_20[names(isoform_sel)[1],], umap_20)
    plot_expression_umaps <- ggplot()+
      geom_point(data = umap_80, aes(x=x,y=y),alpha=0.2,size=0.2, col='grey', show.legend = F)+
      geom_point(data = umap_20, aes(x=x, y=y, col=expr), size=0.7) +
      labs(x="Dim1",y="Dim2",col="scaled expression")+
      scale_colour_gradient2(low = "#313695", mid = "#FFFFBF", high = "#A50026", na.value = NA, midpoint=0)+
      theme_bw()+
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(), 
            panel.border = element_blank(),
            axis.line = element_line(colour = "black"),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            legend.text = element_blank())
    combined_isoform_plot <- plot_grid(plot_isoforms@ggplot,plot_expression_umaps, ncol = 1, rel_heights = c(1, 2.5))
  } else {
    umap_20 <- cbind(t(expr_20), umap_20)
    # ggplot2 lazy evaluation with for loops
    #   may need to unquote idx with !! operator
    plot_idx <- function(idx) {
      p <- ggplot()+
        geom_point(data = umap_80, aes(x=x,y=y),alpha=0.2,size=0.2, col='grey', show.legend = F)+
        geom_point(data = umap_20, aes(x=x,y=y,col=umap_20[,idx]),size=0.7)+
        labs(x="Dim1",y="Dim2",col="scaled expression",title = colnames(umap_20)[idx])+
        scale_colour_gradient2(low = "#313695", mid = "#FFFFBF", high = "#A50026", na.value = NA, midpoint=median(umap_20[,idx]))+
        theme_bw()+
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(), 
              panel.border = element_blank(),
              axis.line = element_line(colour = "black"),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank(),
              legend.position="none")
      p
    }
    plot_expression_umaps <- lapply(1:n_isoforms, plot_idx)
  
    legend <- ggplot(data.frame(x=1:2,y=1:2)) + 
      geom_point(aes(x=x,y=y, col=c(-1,1))) + 
      labs(col="scaled expression") +
      scale_colour_gradient2(low = "#313695", mid = "#FFFFBF", high = "#A50026",
                             breaks = c(-0.8, 0.8),
                             labels=c("Low","High"))
    legend <- get_legend(legend + theme(legend.box.margin = margin(0, 0, 0, 12)))
    combined_isoform_plot <- plot_grid(plotlist =  c(list(annotation=plot_isoforms@ggplot),plot_expression_umaps)) 
    combined_isoform_plot <- plot_grid(combined_isoform_plot, legend, rel_widths = c(3, .4)) 
  }

  if (file.exists(file.path(path,"cluster_annotation.csv")) || (!missing(cluster_annotation) && file.exists(cluster_annotation)) || !is.null(sce_20$cell_type)) {
    cat("Plotting heatmaps ...\n")
    if (file.exists(file.path(path,"cluster_annotation.csv"))) {
      cluster_barcode <- read.csv(file.path(outdir,"cluster_annotation.csv"), stringsAsFactors=FALSE)
      cluster_barcode <- cluster_barcode[,c("barcode_seq", "groups")]
    } else if (!missing(cluster_annotation) && file.exists(cluster_annotation)) {
      cluster_barcode <- read.csv(cluster_annotation, stringsAsFactors=FALSE)
      cluster_barcode <- cluster_barcode[,c("barcode_seq", "groups")]
    } else {
      cluster_barcode <- data.frame(barcode_seq = colnames(sce_20), groups = sce_20$cell_type)
    }
    rownames(cluster_barcode) <- cluster_barcode[,"barcode_seq"] 
    cluster_barcode_20 <- cluster_barcode[rownames(umap_20), "groups", drop=FALSE]
    cluster_barcode_all <- cluster_barcode[colnames(sce_all), "groups", drop=FALSE]
    cell_order <- stats::hclust(stats::dist(t(expr_20)))$order
    heatmap <- pheatmap(expr_20[tr_order,cell_order],
                        cluster_rows = F, cluster_cols = F, 
                        show_colnames = F, show_rownames = T,
                        annotation_col = cluster_barcode_20,
                        width = 4, height = 3)
    heatmap_nolengend <- pheatmap(expr_20[tr_order,cell_order],
                        cluster_rows = F, cluster_cols = F, 
                        show_colnames = F, show_rownames = F,
                        annotation_col = cluster_barcode_20,
                        width = 4, height = 3)
    umap_clusters <- ggplot()+
      geom_point(aes(x=sce_all@int_colData$reducedDims$UMAP[,1], 
                     y=sce_all@int_colData$reducedDims$UMAP[,2], 
                     col=factor(cluster_barcode_all$groups)),
                     size=0.02)+
      labs(x="Dim1", y="Dim2", title = 'Lib_all UMAP', color = "Cluster")
  } else {
    cat("Cluster annotation not found, heatmaps skipped.\n")
    heatmap <- NULL
  }


  outputs <- list(plot_umap = plot_umap,
                plot_isoforms = plot_isoforms, 
                plot_expression_umaps = plot_expression_umaps,
                combined_isoform_plot = combined_isoform_plot) 
  
  if (return_sce_all) {
    outputs[["sce_all"]] <- sce_all
  }

  if (!is.null(heatmap)) {
    outputs[["heatmap"]] <- heatmap
    outputs[["umap_clusters"]] <- umap_clusters
    if (n_isoforms == 2) {
      outputs[["combined_heatmap"]] <- plot_grid(plot_isoforms_plain@ggplot, as.ggplot(heatmap_nolengend), 
                                                 rel_widths = c(1,1.2), align = "hv")
    } else {
      outputs[["combined_heatmap"]] <- plot_grid(plot_isoforms@ggplot, as.ggplot(heatmap_nolengend), 
                                                 rel_widths = c(1,1.2), align = "hv")
    }
  }

  return(outputs)
}
