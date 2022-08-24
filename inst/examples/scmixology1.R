## Process the count matrices from GSE126906 and GSE154869 (scmixology1)
## yields scm_lib80 scm_lib80 scm_lib20

library(SingleCellExperiment)

geodir <- tempfile()
dir.create(geodir)
scmixology1_20_tar <- GEOquery::getGEOSuppFiles(GEO = "GSE154869", baseDir = geodir, fetch_files = TRUE)
utils::untar(rownames(scmixology1_20_tar), exdir = file.path(geodir, "Lib20"))

isoform_gff <- rtracklayer::import.gff3(file.path(geodir, "Lib20", "GSM4681740_isoform_annotated.filtered.gff3.gz"))
isoform_gff$Parent <- as.character(isoform_gff$Parent)
isoform_gff$transcript_id <- unlist(lapply(strsplit(isoform_gff$Parent, split = ":"), function(x) {
    x[2]
}))
isoform_gff <- S4Vectors::split(isoform_gff, isoform_gff$transcript_id)

tr_counts <- read.csv(file.path(geodir, "Lib20", "GSM4681740_transcript_count.csv.gz"))
tr_counts <- tr_counts[tr_counts$transcript_id %in% names(isoform_gff), ]
rownames(tr_counts) <- tr_counts$transcript_id
scm_lib20_transcripts <- SingleCellExperiment(assays = list(counts = tr_counts[, 3:dim(tr_counts)[2]]))
rowRanges(scm_lib20_transcripts) <- isoform_gff[rownames(tr_counts)]
rowData(scm_lib20_transcripts)$gene_id <- gsub("\\..*", "", tr_counts$gene_id)
rowData(scm_lib20_transcripts)$FSM_match <- rownames(scm_lib20_transcripts) # for sc_annotate_plots
rowData(scm_lib20_transcripts)$transcript_id <- rownames(scm_lib20_transcripts)

gene_counts_20 <- read.csv(file.path(geodir, "Lib20", "GSM4681740_gene_count_Lib10.csv.gz"))
rownames(gene_counts_20) <- gene_counts_20$gene_id
gene_counts_20 <- subset(gene_counts_20, select = -c(gene_id))
barcodes_20 <- read.csv(file.path(geodir, "Lib20", "GSM4681740_Lib10.csv.gz"))
colnames(gene_counts_20) <- barcodes_20$barcode_sequence[match(colnames(gene_counts_20), barcodes_20$cell_name)]
scm_lib20 <- SingleCellExperiment(assays = list(counts = gene_counts_20[, 2:dim(gene_counts_20)[2]]))

scmixology1_80_tar <- GEOquery::getGEOSuppFiles(GEO = "GSE126906", baseDir = geodir, fetch_files = TRUE)
utils::untar(rownames(scmixology1_80_tar)[1], exdir = file.path(geodir, "Lib80"))
gene_counts_80 <- read.csv(file.path(geodir, "Lib80", "GSM3618014_gene_count.csv.gz"))
rownames(gene_counts_80) <- gene_counts_80$gene_id
gene_counts_80 <- subset(gene_counts_80, select = -c(gene_id))
barcodes_80 <- read.csv(rownames(scmixology1_80_tar)[2])
colnames(gene_counts_80) <- barcodes_80$barcode_sequence[match(colnames(gene_counts_80), barcodes_80$cell_name)]
scm_lib80 <- SingleCellExperiment(assays = list(counts = gene_counts_80[, 2:dim(gene_counts_80)[2]]))

unlink(geodir, recursive = TRUE)
