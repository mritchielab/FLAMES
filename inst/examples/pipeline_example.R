# download the two fastq files, move them to a folder to be merged together
temp_path <- tempfile()
bfc <- BiocFileCache::BiocFileCache(temp_path, ask=FALSE)
file_url <- 
    "https://raw.githubusercontent.com/OliverVoogd/FLAMESData/master/data"
# download the required fastq files, and move them to new folder
fastq1 <- bfc[[names(BiocFileCache::bfcadd(bfc, "Fastq1", paste(file_url, "fastq/sample1.fastq.gz", sep="/")))]]
fastq2 <- bfc[[names(BiocFileCache::bfcadd(bfc, "Fastq2", paste(file_url, "fastq/sample2.fastq.gz", sep="/")))]]
fastq_dir <- paste(temp_path, "fastq_dir", sep="/") # the downloaded fastq files need to be in a directory to be merged together
dir.create(fastq_dir)
file.copy(c(fastq1, fastq2), fastq_dir)
unlink(c(fastq1, fastq2)) # the original files can be deleted

\dontrun{
# run the FLAMES bulk pipeline, using the downloaded files
se <- bulk_long_pipeline(annot=system.file("extdata/SIRV_anno.gtf", package="FLAMES"), 
                   fastq=fastq_dir,
                   outdir=tempdir(), genome_fa=system.file("extdata/SIRV_genomefa.fasta", package="FLAMES"),
                   config_file=system.file("extdata/SIRV_config_default.json", package="FLAMES"))
}
# OR
# run the FLAMES single cell pipeline
#sce <- sc_long_pipeline(annot, fastq, NULL, outdir, genome_fa, match_barcode=FALSE, config+file=config)