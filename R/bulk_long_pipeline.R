#' Title
#'
#' Semi-supervised isofrom detection and annotation from long read data.
#' This variant is meant for bulk samples. THIS FUNCTION IS GOING TO BE THE ENTIRE 
#' FLAMES PROCESS. IT SHOULD BE BUILT FROM THE SMALLER FUNCTIONS THAT PREVIOUSLY WERE
#' CARRIED OUT UNDER THE HOOD. MOST FUNCTIONS THAT THIS ONE USES SHOULD BE AVAILABLE TO THE USER.
#'
#'  THIS NEED PROPER TYPES
#' @param annot gene annotations in gff3 file format. str
#'
#' @param fastq_dir the folder containing fastq files, each containing data from one sample. str
#'
#' @param in_bam aligned bam file (sorted and indexed). Overwrite --infq? str
#'
#' @param outdir directory to deposite all results in rootdir. Use absolute path? str
#'
#' @param genome_fa genome fasta file. str
#'
#' @param minimap2_dir directory containing minimap2, k8 and paftools.js program. k8 and paftools.js are used to convert gff3 to bed12. str
#'
#' @param config_file json configuration files. str
#'
#' @param downsample_ratio downsampling ratio if performing downsampling analysis. str
#' @importFrom reticulate import_from_path
#' @export
bulk_long_pipeline <- function(annot, fastq_dir, in_bam=NULL, outdir, genome_fa,
                                minimap2_dir=NULL, downsample_ratio,
                                use_junctions=TRUE, no_flank=TRUE,
                                do_genome_align=TRUE, do_isoform_id=TRUE, 
                                do_read_realign=TRUE, do_transcript_quantification=TRUE) {
    # filenames for internal steps
    infq <- paste(outdir, "merged.fastq.gz", sep="/")
    bc_file <- paste(outdir, "pseudo_barcode_annotation.csv", sep="/")

    cat("Preprocessing bulk fastqs...\n")
    # create output directory if one doesn't exist
    if (!dir.exists(outdir)) {
        cat("Output directory does not exists: one is being created\n")
        dir.create(outdir)
        print(outdir)
    }
    # run the merge_bulk_fastq function as preprocessing
    merge_bulk_fastq(fastq_dir, bc_file, infq)

    cat("Running FLAMES pipeline...\n")
    # argument verificiation
    if (downsample_ratio > 1 || downsample_ratio <= 0) {
        stop("downsample_ratio should be between 0 and 1")
    }
    if (!file.exists(infq) || !file.exists(annot) || !file.exists(genome_fa)) {
        stop(paste0("Make sure all files exists: ", infq, ", ", annot, ", ", genome_fa))
    }
    if (is.null(in_bam)) {
        in_bam = ""
    } else {
        if (!file.exists(in_bam)) {
            stop("Make sure input in_bam file exists")
        }
    }
    if (is.null(minimap2_dir)) minimap2_dir = ""

    # setup of internal arguments which hold output files and intermediate files
    isoform_gff3 = paste(outdir, "isoform_annotated.gff3", sep="/")
    isoform_gff3_f = paste(outdir, "isoform_annotated.filtered.gff3", sep="/")
    FSM_anno_out = paste(outdir, "isoform_FSM_annotation.csv", sep="/")
    raw_splice_isoform = paste(outdir, "splice_raw.gff3", sep="/")
    tss_tes_stat = paste(outdir, "tss_tes.bedgraph", sep="/")
    transcript_fa = paste(outdir, "transcript_assembly.fa", sep="/")
    transcript_fa_idx = paste(outdir, "transcript_assembly.fa.fai", sep="/")
    tmp_bam = paste(outdir, "tmp.align.bam", sep="/")
    tmp_bed = paste(outdir, "tmp.splice_anno.bed12", sep="/")
    genome_bam = paste(outdir, "align2genome.bam", sep="/")
    realign_bam = paste(outdir, "realign2transcript.bam", sep="/")
    tr_cnt_csv = paste(outdir, "transcript_count.csv.gz", sep="/")
    tr_badcov_cnt_csv = paste(outdir, "transcript_count.bad_coverage.csv.gz", sep="/")

    cat("Input parameters:\n")
    cat("\tgene annotation:", annot, "\n")
    cat("\tgenome fasta:", genome_fa, "\n")
    if (in_bam != "") {
        cat("\tinput bam:", in_bam, "\n")
        genome_bam = in_bam
    } else cat("\tinput fastq:", infq, "\n")
    cat("\toutput directory:", outdir, "\n")
    cat("\tdirectory containing minimap2:", minimap2_dir, "\n")

    # align reads to genome
    if (in_bam=="" && do_genome_align) {
        cat("Aligning reads to genome using minimap2\n")
        align_reads_to_genome(annot, infq, genome_fa, genome_bam, minimap2_dir, outdir, use_junctions, no_flank)
    } else {
        cat("Skip aligning reads to genome\n")
    }

    
    # find isofrom


    invisible()
}