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
#' @param path The path to the output folder of \code{sc_long_pipeline}
#' the folder needs to contain:
#' \itemize{
#'  \item{transcript_count.csv.gz}{ - the transcript count matrix }
#'  \item{isoform_FSM_annotation.csv}{ - the full splice match annotation file }
#'  \item{cluster_annotation.csv}{ - cluster annotation file }
#' }
#'
#' @param min_count The minimum UMI count threshold for filtering isoforms.
#'
#' @return a \code{data.frame} containing the following columns:
#' \itemize{
#'  \item{gene_name}{ - differentially transcribed genes }
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
#' @importFrom dplyr group_by summarise_at top_n left_join summarise groups
#' @importFrom tidyr gather pivot_wider as_tibble
#' @importFrom magrittr "%>%"
#' @importFrom S4Vectors DataFrame
#' @importFrom SingleCellExperiment counts SingleCellExperiment
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
#'     reference_csv = bc_allow
#'   )
#'   group_anno <- data.frame(barcode_seq = colnames(sce), groups = SingleCellExperiment::counts(sce)["ENSMUST00000169826.2", ] > 1)
#'   write.csv(group_anno, file.path(outdir, "cluster_annotation.csv"), row.names = FALSE)
#'   sc_DTU_analysis(sce, outdir, min_count = 1)
#' }
sc_DTU_analysis <- function(sce, path, min_count = 15) {

  # Check input

  if (!missing(sce) && !missing(path)) {
    cat("Both sce and path arguments are provided, ignoring the path argument.\n")
  }

  # sce object from sc_long_pipeline
  if (!missing(sce)) {
    if (!is(sce, "SingleCellExperiment")) {
      stop("sce need to be an SingleCellExperiment Object returned by sc_long_pipeline()")
    }
    if (!file.exists(file.path(sce@metadata$OutputFiles$outdir, "isoform_FSM_annotation.csv"))) {
      stop("Missing isoform_FSM_annotation.csv")
    }
    if (!file.exists(file.path(sce@metadata$OutputFiles$outdir, "cluster_annotation.csv"))) {
      stop(paste(c("Please provide a cluster_annotation.csv under ", sce@metadata$OutputFiles$outdir)))
    }

    outdir <- sce@metadata$OutputFiles$outdir
    transcript_count <- cbind(rowData(sce), as.data.frame(counts(sce)))

    if (!missing(path) && (path != outdir)) {
      cat("Warning: the output files of the SCE object are in")
      cat(sce@metadata$OutputFiles$outdir)
      cat(", which is different from the path argument.\n")
    }
  } else if (is.character(path)) {
    # Path to output dir

    if (!file.exists(file.path(path, "transcript_count.csv.gz"))) {
      stop(paste(c("transcript_count.csv.gz not found under ", path)))
    }
    if (!file.exists(file.path(path, "isoform_FSM_annotation.csv"))) {
      stop(paste(c("isoform_FSM_annotation.csv not found under ", path)))
    }
    if (!file.exists(file.path(path, "cluster_annotation.csv"))) {
      stop(paste(c("cluster_annotation.csv not found under ", path)))
    }

    outdir <- path
    cat("Loading transcript_count.csv.gz ...\n")
    transcript_count <- read.csv(file.path(outdir, "transcript_count.csv.gz"), stringsAsFactors = FALSE)
  } else {
    stop("Please provide the sce arugment or the path argument")
  }

  cat("Loading isoform_FSM_annotation.csv ...\n\n")
  isoform_FSM_annotation <- read.csv(file.path(outdir, "isoform_FSM_annotation.csv"), stringsAsFactors = FALSE)
  cat("Selecting transcript_ids with full splice match...\n")
  transcript_count <- transcript_count[match(isoform_FSM_annotation$transcript_id, transcript_count$transcript_id), ]
  transcript_count$FSM_match <- isoform_FSM_annotation$FSM_match
  cat("Summing transcripts with same FSM ... \n")
  cell_bcs <- colnames(transcript_count)[!(colnames(transcript_count) %in% c("transcript_id", "gene_id", "FSM_match"))]

  # mer_tmp: cell_bcs as cols, FSM_match as rows.
  fsm_csv <- file.path(outdir, "FSM_count.csv.gz")
  if (file.exists(fsm_csv)) {
    mer_tmp <- read.csv(fsm_csv)
  } else {
    # sum transcript (FSM) counts
    mer_tmp <- as_tibble(transcript_count) %>%
      group_by(FSM_match) %>%
      summarise_at(cell_bcs, sum)
    cat("Creating FSM_count.csv.gz ...\n")
    write.csv(mer_tmp, file = gzfile(fsm_csv), row.names = FALSE)
    cat(paste0(c(fsm_csv, "saved.\n")))
  }

  tr_anno <- transcript_count[, c("transcript_id", "gene_id", "FSM_match")]
  rm(transcript_count)
  # keep only one transcript_id for each FSM_match
  tr_anno <- tr_anno[match(mer_tmp$FSM_match, tr_anno$FSM_match), ]

  cat("QC & outlier removal: \n")
  tr_sce <- SingleCellExperiment(assays = list(counts = as.matrix(mer_tmp[, -1])))
  rownames(tr_sce) <- mer_tmp$FSM_match
  rowData(tr_sce) <- DataFrame(tr_anno)
  tr_sce <- addPerCellQC(tr_sce) # Consider allowing user to change parameters for QC?
  tr_sce <- addPerFeatureQC(tr_sce)
  keep.hi <- isOutlier(tr_sce$sum, type = "higher", log = TRUE)
  keep.low <- isOutlier(tr_sce$sum, type = "lower", log = TRUE)
  tr_sce <- tr_sce[, (!keep.hi) & (!keep.low)]
  cat(paste(c(sum(!((!keep.hi) & (!keep.low))), " cell BC(s) removed, ", dim(tr_sce)[2], " remaining\n")))

  # Remove version number from gene_id
  rowData(tr_sce)$gene_id <- gsub("\\..*", "", rowData(tr_sce)$gene_id)
  # gene_name = mapIds(org.Hs.eg.db,
  #                   keys=rowData(tr_sce)$gene_id,
  #                   column="SYMBOL",
  #                   keytype="ENSEMBL",
  #                   multiVals="first")

  # gene_name[is.na(gene_name)] = rowData(tr_sce)$gene_id[is.na(gene_name)]
  rowData(tr_sce)$gene_name <- rowData(tr_sce)$gene_id

  cat("Loading cluster_annotation.csv ...\n")
  cluster_barcode_anno <- read.csv(file.path(outdir, "cluster_annotation.csv"), stringsAsFactors = FALSE)
  rownames(cluster_barcode_anno) <- cluster_barcode_anno$barcode_seq
  comm_cells <- intersect(colnames(tr_sce), rownames(cluster_barcode_anno))
  tr_sce <- tr_sce[, comm_cells]

  cat("Adding cluster information to SCE object ...\n")
  cluster_barcode_anno <- cluster_barcode_anno[comm_cells, ]
  colData(tr_sce) <- cbind(colData(tr_sce), DataFrame(cluster_barcode_anno))

  cat("Filtering for genes with at least 2 detected isforms ...")
  tr_sce_multi <- tr_sce[rowData(tr_sce)$gene_name %in% names(table(rowData(tr_sce)$gene_name)[table(rowData(tr_sce)$gene_name) > 1]), ]
  cat(paste(c("     ", dim(tr_sce_multi)[1], "FSM_match(s) left.\n")))

  cat("Keeping only the top 4 expressed FSM_matches for each gene ...")
  row_meta <- as.data.frame(rowData(tr_sce_multi))
  row_meta <- row_meta %>%
    group_by(gene_name) %>%
    top_n(n = 4, wt = mean) # Consider adding a parameter here to allow user to decide how many to keep?
  cat(paste(c("     ", dim(row_meta)[1], "FSM_match(s) left.\n")))

  # Apply the top_n filtering to the sce object
  tr_sce_multi <- tr_sce_multi[rowData(tr_sce_multi)$transcript_id %in% row_meta$transcript_id, ]

  tmp_df <- as.data.frame(counts(tr_sce_multi))
  tmp_df$tr_id <- rownames(tmp_df)
  tr_sce_multi$barcode_seq <- colnames(tr_sce_multi)
  data_long <- gather(tmp_df, cell_id, cnt, colnames(tmp_df)[!(colnames(tmp_df) == "tr_id")])
  tr_sce_multi$groups <- as.factor(tr_sce_multi$groups)
  data_long <- left_join(data_long, as.data.frame(colData(tr_sce_multi)[, c("barcode_seq", "groups")]), by = c("cell_id" = "barcode_seq"))
  data_long <- data_long %>%
    group_by(tr_id, groups) %>%
    summarise(cnt = sum(cnt))
  data_wide <- data_long %>% pivot_wider(id_cols = tr_id, names_from = groups, values_from = cnt)
  all_grp <- colnames(data_wide)[-1]
  # data_wide = as.matrix(data_wide[,-1])
  row_meta <- as.data.frame(rowData(tr_sce_multi))
  row_meta <- row_meta %>% left_join(data_wide, by = c("FSM_match" = "tr_id"))

  # Filter for genes of which the most aboundant isoform have at least 10 counts in one cell group (?)
  # Keep only 2 most aboundant isoforms of each gene for each cell group (?)
  cat("Filtering isoforms ... ")
  get_rm <- function(x) {
    tmp <- row_meta
    tmp$l <- tmp[, x]
    gc <- tmp %>%
      group_by(gene_name) %>%
      dplyr::summarise(l_s = max(l))
    tmp <- tmp[tmp$gene_name %in% gc$gene_name[gc$l_s > 10], ]
    tmp <- tmp %>%
      group_by(gene_name) %>%
      top_n(n = 2, wt = l)
    return(tmp$transcript_id)
  }
  sel_tr <- Reduce(union, lapply(all_grp, get_rm))
  row_meta <- row_meta[row_meta$transcript_id %in% sel_tr, ]
  row_meta <- row_meta %>%
    group_by(gene_name) %>%
    top_n(n = 4, wt = mean)
  cat(paste(c(dim(row_meta)[1], " transcript_id(s) remaining.\n")))

  cat("Performing Chi-square tests ...\n")
  ge_name <- c()
  X_value <- c()
  df <- c()
  p_value <- c()
  DTU_tr <- c()
  DTU_group <- c()
  for (ge in unique(row_meta$gene_name)) {
    # for (ge in sel_ge){
    data_wide_tmp <- data_wide[rowData(tr_sce_multi)$gene_name == ge, ]
    data_wide_tmp <- data_wide_tmp[data_wide_tmp$tr_id %in% row_meta$FSM_match, ]
    rn <- data_wide_tmp$tr_id
    data_wide_tmp <- as.matrix(data_wide_tmp[, -1])
    rownames(data_wide_tmp) <- rn
    data_wide_tmp <- data_wide_tmp[apply(data_wide_tmp, 1, max) > min_count, apply(data_wide_tmp, 2, max) > min_count]
    if (!is.null(dim(data_wide_tmp))) {
      if (ncol(data_wide_tmp) >= 2 & nrow(data_wide_tmp) > 1) {
        fit <- suppressWarnings(chisq.test(data_wide_tmp))
        ge_name <- c(ge_name, ge)
        X_value <- c(X_value, fit$statistic)
        df <- c(df, fit$parameter)
        p_value <- c(p_value, fit$p.value)
        cs <- colSums(fit$residuals^2)
        DTU_group <- c(DTU_group, names(cs[order(cs, decreasing = TRUE)[1]]))
        if (nrow(data_wide_tmp) == 2) {
          rs <- rowSums(data_wide_tmp)
          DTU_tr <- c(DTU_tr, names(rs[order(rs)[1]]))
        } else {
          hi1 <- which(rowSums(data_wide_tmp) == max(rowSums(data_wide_tmp)))
          if (length(hi1) > 1) {
            re1 <- rowSums(fit$residuals[hi1, ]^2)
          } else {
            re1 <- rowSums(fit$residuals[-hi1, ]^2)
          }
          DTU_tr <- c(DTU_tr, names(re1[order(re1, decreasing = TRUE)[1]]))
        }
      }
    }
  }
  res_df <- data.frame(
    gene_name = ge_name,
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
  write.csv(res_df, file = file.path(path, "sc_DTU_analysis.csv"), row.names = FALSE)
  cat(paste(c("Results saved to ", file.path(path, "sc_DTU_analysis.csv"), "\n")))
  return(res_df)
}
