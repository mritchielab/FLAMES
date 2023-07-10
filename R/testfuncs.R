#' @export
testfuncR <- function() {
    # data <- "/stornext/General/data/user_managed/grpu_mritchie_1/Oliver/FLAMESTesting/FLAMESData/FLAMESbulkdata/SIRV/"
    data <- "/Volumes/MattLab/Oliver/FLAMESTesting/FLAMESData/FLAMESbulkdata/SIRV/"
    genome_bam <- paste0(data, "FLAMESout/align2genome.bam")
    gff3 <- paste0(data,"SIRV_isoforms_multi-fasta-annotation_C_170612a.gtf")
    genomefa <- paste0(data, "SIRV_isoforms_multi-fasta_170612a.fasta")
    config_file <- paste0(data, "SIRV_config_default.json")
    
    out = paste0(tempdir(), "/")
    isoform_gff3 <- paste0(out, "isoform_annotated.gff3")
    tss_tes_stat <- paste0(out, "tss_tes.bedgraph")
    transcript_fa <- paste0(out, "transcript_assembly.fa")

    config <- jsonlite::fromJSON(config_file)

    file.exists(c(genome_bam, gff3, genomefa, config_file))
    testfunc(gff3, genome_bam, isoform_gff3, tss_tes_stat, genomefa, transcript_fa, config$isoform_parameters, "")

    print(out)
}