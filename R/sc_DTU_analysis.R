# TODO:
#       Add parameters
#
#
#' FLAMES Differential Transcript Usage Analysis
#'
#' Chi-square based differential transcription usage analysis. This variant is meant for single cell data.
#' Takes the \code{SingleCellExperiment} object from \code{sc_long_pipeline} as input.
#' Alternatively, the path to the output folder could be provided instead of the SCE object.
#' A cluster annotation file \code{cluster_annotation.csv} is required, please provide this file under the
#' output folder of \code{sc_long_pipeline}.
#'
#'
#' @details
#' This function will search for genes that have at least two isoforms, each with more than \code{min_count} UMI counts.
#' For each gene, the per cell transcript counts were merged by group to generate pseudo bulk samples.
#' Grouping is specified by the \code{cluster_annotation.csv} file.
#' The top 2 highly expressed transcripts for each group were selected and a UMI count matrix where
#' the rows are selected transcripts and columns are groups was used as input to a chi-square test of independence (chisq.test).
#' Adjusted P-values were calculated by Benjamini–Hochberg correction.
#'
#' @param sce The \code{SingleCellExperiment} object from \code{sc_long_pipeline}, an additional \code{cluster_annotation.csv}
#' file is required under the output folder of the SCE object.
#'
#' @param min_count The minimum UMI count threshold for filtering isoforms.
#'
#' @return a \code{data.frame} containing the following columns:
#' \itemize{
#'  \item{gene_id}{ - differentially transcribed genes }
#'  \item{X_value}{ - the X value for the DTU gene}
#'  \item{df}{ - degrees of freedom of the approximate chi-squared distribution of the test statistic }
#'  \item{DTU_tr}{ - the transcript_id with the highest squared residuals}
#'  \item{DTU_group}{ - the cell group with the highest squared residuals}
#'  \item{p_value}{ - the p-value for the test}
#'  \item{adj_p}{ - the adjusted p-value (by Benjamini–Hochberg correction)}
#' }
#' The table is sorted by decreasing P-values. It will also be saved as \code{sc_DTU_analysis.csv} under the
#' output folder.
#'
#' @importFrom dplyr group_by summarise_at top_n left_join summarise groups mutate filter_at any_vars
#' @importFrom tidyr gather pivot_wider as_tibble
#' @importFrom magrittr "%>%"
#' @importFrom S4Vectors DataFrame
#' @importFrom SingleCellExperiment counts SingleCellExperiment colLabels colLabels<-
#' @importFrom SummarizedExperiment rowData colData rowData<- colData<-
#' @importFrom scuttle addPerCellQC addPerFeatureQC isOutlier
#' @importFrom utils write.csv
#' @importFrom stats chisq.test
#' @export
#' @examples
#' outdir <- tempfile()
#' dir.create(outdir)
#' bc_allow <- file.path(outdir, "bc_allow.tsv")
#' genome_fa <- file.path(outdir, "rps24.fa")
#' R.utils::gunzip(filename = system.file("extdata/bc_allow.tsv.gz", package = "FLAMES"), destname = bc_allow, remove = FALSE)
#' R.utils::gunzip(filename = system.file("extdata/rps24.fa.gz", package = "FLAMES"), destname = genome_fa, remove = FALSE)
#'
#' if (is.character(locate_minimap2_dir())) {
#'   sce <- FLAMES::sc_long_pipeline(
#'     genome_fa = genome_fa,
#'     fastq = system.file("extdata/fastq", package = "FLAMES"),
#'     annotation = system.file("extdata/rps24.gtf.gz", package = "FLAMES"),
#'     outdir = outdir,
#'     match_barcode = TRUE,
#'     barcodes_file = bc_allow
#'   )
#'   group_anno <- data.frame(barcode_seq = colnames(sce), groups = SingleCellExperiment::counts(sce)["ENSMUST00000169826.2", ] > 1)
#'   write.csv(group_anno, file.path(outdir, "cluster_annotation.csv"), row.names = FALSE)
#'   sc_DTU_analysis(sce, min_count = 1)
#' }
sc_DTU_analysis <- function(sce, min_count = 15) {

  # sce object from sc_long_pipeline
  if (!is(sce, "SingleCellExperiment")) {
    stop("sce need to be an SingleCellExperiment Object returned by sc_long_pipeline()")
  }
  if (!file.exists(file.path(sce@metadata$OutputFiles$outdir, "isoform_FSM_annotation.csv"))) {
    stop("Missing isoform_FSM_annotation.csv")
  }
  stopifnot(!is.null(colLabels(sce)) || 
              file.exists(file.path(sce@metadata$OutputFiles$outdir, "cluster_annotation.csv")))

  outdir <- sce@metadata$OutputFiles$outdir
  cat("Loading isoform_FSM_annotation.csv ...\n")
  isoform_FSM_annotation <- read.csv(file.path(outdir, "isoform_FSM_annotation.csv"), stringsAsFactors = FALSE)
  
  cat("Selecting transcript_ids with full splice match ...\n")
  sce <- sce[rownames(sce) %in% isoform_FSM_annotation$transcript_id,]
  rowData(sce)$FSM_match <- 
    isoform_FSM_annotation[
      match(rownames(sce),isoform_FSM_annotation$transcript_id),
      "FSM_match"
      ]

  cat("Summing transcripts with same FSM ... \n")
  # mer_tmp: cell_bcs as cols, FSM_match as rows.
  fsm_csv <- file.path(outdir, "FSM_count.csv.gz")
  if (file.exists(fsm_csv)) {
    mer_tmp <- read.csv(fsm_csv)
  } else {
    # sum transcript (FSM) counts
    mer_tmp <- as_tibble(counts(sce)) %>%
      mutate(FSM_match = rowData(sce)$FSM_match) %>%
      group_by(FSM_match) %>%
      summarise_at(colnames(sce), sum)
    cat("Creating FSM_count.csv.gz ...\n")
    write.csv(mer_tmp, file = gzfile(fsm_csv), row.names = FALSE)
    cat(paste0(c(fsm_csv, "saved.\n")))
  }
  stopifnot(all(complete.cases(mer_tmp)))
  cat("\t", nrow(mer_tmp), "FSM_match(s) found.\n")

  tr_sce <- SingleCellExperiment(assays = list(counts = as.matrix(mer_tmp[, -1])))
  rownames(tr_sce) <- mer_tmp$FSM_match
  rowData(tr_sce) <- rowData(sce)[match(mer_tmp$FSM_match, rowData(sce)$FSM_match),]

  # Remove version number from gene_id
  rowData(tr_sce)$gene_id <- gsub("\\..*", "", rowData(tr_sce)$gene_id)

  if (is.null(colLabels(sce))) {
    cat("Loading cluster_annotation.csv ...\n")
    cluster_barcode_anno <- read.csv(file.path(outdir, "cluster_annotation.csv"), stringsAsFactors = FALSE)
    comm_cells <- na.omit(match(colnames(tr_sce), cluster_barcode_anno$barcode_seq))
    cluster_barcode_anno <- cluster_barcode_anno[comm_cells,]
    
    cat("Adding cluster information to SCE object ...\n")
    colLabels(tr_sce) <- cluster_barcode_anno$group
  } else {
    colLabels(tr_sce) <- colLabels(sce)
  }

  cat("Filtering for genes with at least 2 detected isforms ...")
  tr_sce_multi <- tr_sce[rowData(tr_sce)$gene_id %in% names(table(rowData(tr_sce)$gene_id)[table(rowData(tr_sce)$gene_id) > 1]), ]
  cat(paste(c("     ", dim(tr_sce_multi)[1], "FSM_match(s) left.\n")))

  cat("Keeping only the top 4 expressed FSM_matches for each gene ...")
  rowData(tr_sce_multi)$mean <- rowMeans(counts(tr_sce_multi))
  top_4s <- rowData(tr_sce_multi) %>%
    as.data.frame() %>%
    group_by(gene_id) %>%
    top_n(n = 4, wt = mean) # Consider adding a parameter here to allow user to decide how many to keep?
  cat(paste(c("     ", dim(top_4s)[1], "FSM_match(s) left.\n")))

  # Apply the top_n filtering to the sce object
  tr_sce_multi <- tr_sce_multi[rowData(tr_sce_multi)$transcript_id %in% top_4s$transcript_id, ]
  
  cat("Aggregating counts by cluster labels ...\n")
  tr_sce_multi$barcode_seq <- colnames(tr_sce_multi)
  counts_group <- counts(tr_sce_multi) %>%
    as.data.frame() %>%
    mutate(tr_id = rownames(.)) %>%
    gather(cell_id, cnt, -"tr_id") %>% # long format: tr_id, cell_id, cnt
    left_join(as.data.frame(colData(tr_sce_multi)[, c("barcode_seq", "label")]),
              by = c("cell_id" = "barcode_seq")) %>%
    group_by(tr_id, label) %>%
    summarise(cnt = sum(cnt)) %>%
    pivot_wider(id_cols = tr_id, names_from = label, values_from = cnt) %>%
    left_join(as.data.frame(rowData(tr_sce_multi)), by = c("tr_id" = "FSM_match"))
  stopifnot(all(complete.cases(counts_group)))
  # counts_group: transcript counts by group
  # tr_id, cluster_1, cluster2, ...

  cat("Filtering isoforms ... ")
  filter_tr <- function(x) {
    # for cluster x, filter genes for at least min_count for top expressed isoform
    group_max <- counts_group %>%
      group_by(gene_id) %>%
      summarise_at(x, max) %>%
      filter_at(x, all_vars(. > min_count)) 
    
    # select the top 2 isoforms expressed by group x for each gene
    filtered_tr_ids <- counts_group %>% 
      filter_at("gene_id", all_vars(. %in% group_max$gene_id)) %>%
      group_by(gene_id) %>%
      top_n(n = 2, wt = x) %>%
      dplyr::pull("tr_id")
    return(filtered_tr_ids)
  }
  
  all_grp <- levels(colLabels(tr_sce_multi))
  sel_tr <- Reduce(union, lapply(all_grp, filter_tr))
  counts_group <- counts_group[counts_group$tr_id %in% sel_tr, ]
  counts_group <- counts_group %>%
    group_by(gene_id) %>%
    top_n(n = 4, wt = mean)
  cat(paste(c(nrow(counts_group), " transcript_id(s) remaining.\n")))

  cat("Performing Chi-square tests ...\n")
  ge_name <- c()
  X_value <- c()
  df <- c()
  p_value <- c()
  DTU_tr <- c()
  DTU_group <- c()
  pb_max <- length(unique(counts_group$gene_id))
  if (pb_max > 10) {
    pb <- txtProgressBar(min = 1, max = length(unique(counts_group$gene_id)),
                        initial = 1, style=3)
    pbi <- 0
  }
  for (gene in unique(counts_group$gene_id)) {
    # transcript x clutser count matrix for gene
    counts_group_gene <- counts_group %>%
      filter_at("gene_id", all_vars(. == gene)) %>%
      filter_at(all_grp, any_vars(. > min_count))
    mtx <- counts_group_gene %>%
      ungroup() %>%
      select(all_of(all_grp)) %>%
      select_if(function(col_x){max(col_x) > min_count}) %>%
      as.matrix()
    rownames(mtx) <- counts_group_gene$tr_id
    
    if (!is.null(dim(mtx))) {
      if (ncol(mtx) > 1 && nrow(mtx) > 1) {
        fit <- suppressWarnings(chisq.test(mtx))
        ge_name <- c(ge_name, gene)
        X_value <- c(X_value, fit$statistic)
        df <- c(df, fit$parameter)
        p_value <- c(p_value, fit$p.value)
        cs <- colSums(fit$residuals^2)
        DTU_group <- c(DTU_group, names(cs[order(cs, decreasing = TRUE)[1]]))
        if (nrow(mtx) == 2) {
          rs <- rowSums(mtx)
          DTU_tr <- c(DTU_tr, names(rs[order(rs)[1]]))
        } else {
          hi1 <- which(rowSums(mtx) == max(rowSums(mtx)))
          if (length(hi1) > 1) {
            re1 <- rowSums(fit$residuals[hi1, ]^2)
          } else {
            re1 <- rowSums(fit$residuals[-hi1, ]^2)
          }
          DTU_tr <- c(DTU_tr, names(re1[order(re1, decreasing = TRUE)[1]]))
        }
      }
    }
    if (exists("pb")) {
      pbi <- pbi + 1
      setTxtProgressBar(pb,pbi)
    }
  }
  if (exists("pb")) {
    close(pb)
  }
  res_df <- data.frame(
    gene_id = ge_name,
    X_value = X_value,
    df = df,
    DTU_tr = DTU_tr,
    DTU_group = DTU_group,
    p_value = p_value, stringsAsFactors = FALSE
  )
  if (any(dim(res_df) == 0)) {
    message("No DTU gene was found\n")
    return(res_df)
  }
  res_df$adj_p <- res_df$p_value * nrow(res_df)
  res_df$adj_p <- sapply(res_df$adj_p, function(x) {
    min(1, x)
  })
  res_df <- res_df[order(res_df$p_value), ]
  warning("Chi-squared approximation(s) may be incorrect")
  write.csv(res_df, file = file.path(outdir, "sc_DTU_analysis.csv"), row.names = FALSE)
  cat(paste(c("Results saved to ", file.path(outdir, "sc_DTU_analysis.csv"), "\n")))
  return(res_df)
}
