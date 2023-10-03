## Process the count matrices from GSE126906 and GSE154869
## yields scmixology_lib90 scmixology_lib10 scmixology_lib10_transcripts

library(SingleCellExperiment)

geodir <- tempfile()
dir.create(geodir)
scmixology1_90_tar <- GEOquery::getGEOSuppFiles(GEO = "GSE154869", baseDir = geodir, fetch_files = TRUE)
utils::untar(rownames(scmixology1_90_tar), exdir = file.path(geodir, "Lib90"))

# scmixology_lib10_transcripts
isoform_gff <- rtracklayer::import.gff3(file.path(geodir, "Lib90", "GSM4681740_isoform_annotated.filtered.gff3.gz"))
isoform_gff$Parent <- as.character(isoform_gff$Parent)
isoform_gff$transcript_id <- unlist(lapply(strsplit(isoform_gff$Parent, split = ":"), function(x) {
    x[2]
}))
isoform_gff <- S4Vectors::split(isoform_gff, isoform_gff$transcript_id)

tr_counts <- read.csv(file.path(geodir, "Lib90", "GSM4681740_transcript_count.csv.gz"))
tr_counts <- tr_counts[tr_counts$transcript_id %in% names(isoform_gff), ]
rownames(tr_counts) <- tr_counts$transcript_id
scmixology_lib10_transcripts <- SingleCellExperiment(assays = list(counts = tr_counts[, 3:dim(tr_counts)[2]]))
rowRanges(scmixology_lib10_transcripts) <- isoform_gff[rownames(tr_counts)]
rowData(scmixology_lib10_transcripts)$gene_id <- gsub("\\..*", "", tr_counts$gene_id)
rowData(scmixology_lib10_transcripts)$FSM_match <- rownames(scmixology_lib10_transcripts) # for sc_annotate_plots
rowData(scmixology_lib10_transcripts)$transcript_id <- rownames(scmixology_lib10_transcripts)

# scmixology_lib10
gene_counts_90 <- read.csv(file.path(geodir, "Lib90", "GSM4681740_gene_count_Lib10.csv.gz"))
rownames(gene_counts_90) <- gene_counts_90$gene_id
gene_counts_90 <- subset(gene_counts_90, select = -c(gene_id))
barcodes_90 <- read.csv(file.path(geodir, "Lib90", "GSM4681740_Lib10.csv.gz"))
colnames(gene_counts_90) <- barcodes_90$barcode_sequence[match(colnames(gene_counts_90), barcodes_90$cell_name)]
scmixology_lib10 <- SingleCellExperiment(assays = list(counts = gene_counts_90[, 2:dim(gene_counts_90)[2]]))

# scmixology_lib90
scmixology1_90_tar <- GEOquery::getGEOSuppFiles(GEO = "GSE126906", baseDir = geodir, fetch_files = TRUE)
utils::untar(rownames(scmixology1_90_tar)[1], exdir = file.path(geodir, "Lib90"))
gene_counts_90 <- read.csv(file.path(geodir, "Lib90", "GSM3618014_gene_count.csv.gz"))
rownames(gene_counts_90) <- gene_counts_90$gene_id
gene_counts_90 <- subset(gene_counts_90, select = -c(gene_id))
barcodes_90 <- read.csv(rownames(scmixology1_90_tar)[2])
colnames(gene_counts_90) <- barcodes_90$barcode_sequence[match(colnames(gene_counts_90), barcodes_90$cell_name)]
scmixology_lib90 <- SingleCellExperiment(assays = list(counts = gene_counts_90[, 2:dim(gene_counts_90)[2]]))

unlink(geodir, recursive = TRUE)

# filter rows and cols to reduce size
scmixology_lib90 <- scuttle::addPerCellQC(scmixology_lib90)
scmixology_lib90 <- scuttle::addPerFeatureQC(scmixology_lib90)
# scmixology_lib10 <- scuttle::addPerCellQC(scmixology_lib10)
# scmixology_lib10 <- scuttle::addPerFeatureQC(scmixology_lib10)
scmixology_lib90 <- scmixology_lib90[(rowData(scmixology_lib90)$detected > 5 & rowData(scmixology_lib90)$mean > 2),]

# WIP

#gene_intersection <- intersect(rownames(scmixology_lib90), rownames(scmixology_lib10))
#scmixology_lib90 <- scmixology_lib90[gene_intersection, ]
#scmixology_lib10 <- scmixology_lib10[gene_intersection, ]
