#' GFF3 to BED12
#'
#' Converts a gff3 file to a bed12, using paftools.js from the minimap2 package
#' 
#' @param minimap2_prog_path Absolute path to the directory containing k8, paftools,js and minimap2 executables from the minimap2 package. Leave as default if these are in the current working directory
#' @param gff3_file The gff3_file to convert
#' @param bed12_file The filename of the bed12 output file.
#' 
#' @return file path to the created bed12_file
#' @importFrom reticulate import_from_path
#' 
#' @export
gff3_to_bed12 <- function(minimap2_prog_path=NULL, gff3_file, bed12_file) {
    python_path <- system.file("python", package="FlamesR")
    callBasilisk(flames_env, function(mm2_path, gff3, bed12) {
        align <-reticulate::import_from_path("minimap2_align", python_path)

        if (is.null(mm2_path)) mm2_path = ""
        align$gff3_to_bed12(mm2_path, gff3, bed12)
    }, mm2_path=minimap2_prog_path, gff3=gff3_file, bed12=bed12_file)

    bed12_file # output file
}

#' Minimap2 Align to Genome
#'
#' @description
#' Uses minimap2 to align sequences agains a reference databse. 
#' Uses options "-ax splice -t 12 -k14 --secondary=no \code{fa_file} \code{fq_in} | samtools view -bS -@ 4 -m 2G -o \code{bam_out}"
#' 
#' @param minimap2_prog_path Absolute path to the directory containing minimap2
#' @param fa_file Fasta file used as a reference database for alignment
#' @param fq_in Fastq file used as a query sequence file
#' @param bam_out Output BAM file
#' @param no_flank Boolean; used if studying SIRV, to let minimap2 ignore additional bases
#' @param bed12_junc Gene annotations in BED12 format. If specified, minmap2 prefers splicing in annotations.
#'
#' @return file path to the given output BAM file, \code{bam_out}
#' @importFrom reticulate import_from_path
#'
#' @export
minimap2_align <- function(minimap2_prog_path=NULL, fa_file, fq_in, bam_out, no_flank=FALSE, bed12_junc=NULL) {
    callBasilisk(flames_env, function (mm2_path, fa, fq, bam, flank, bed12_junc) {
        python_path <- system.file("python", package="FlamesR")
        mm2 <- reticulate::import_from_path("minimap2_align", python_path)

        if (is.null(minimap2_prog_path)) minimap2_prog_path = ""
        mm2$minimap2_align(mm2_path, fa, fq, bam, flank, bed12_junc)
    }, mm2_path=minimap2_prog_path, fa=fa_file, fq=fq_in, bam=bam_out, flank=no_flank, bed12_junc=bed12_junc)
    bam_out # output file
}

#' Samtools Sort Index
#'
#' Sorts and then indexes alignments for fast random access using samtools. Uses sortBam and indexBam from
#' Rsamtools.
#' 
#' @param bam_in Input BAM File
#' @param bam_out Output BAM file
#' @importFrom Rsamtools sortBam indexBam
#' @return the path to the output file, given as \code{bam_out}
#' @export
samtools_sort_index <- function(bam_in, bam_out) {
    sortBam(bam_in, gsub("\\.bam", "", bam_out))
    indexBam(bam_out)

    bam_out
}

#' Minimap2 Align to Transcript
#'
#' @description
#' Uses minimap2 to align to transcript. 
#' Uses options "-ax map-ont -p 0.9 --end-bonus 10 -N 3 -t 12 \code{fa_file} \code{fq_in} | samtools view -bS -@ 4 -m 2G -o \code{bam_out}"
#' 
#' @param mm2_prog_path Absolute path to the directory containing minimap2
#' @param fa_file Input fasta file used as a reference database
#' @param fq_in Input fastq used as a query sequence file
#' @param bam_out Output BAM file, containing aligned sequences
#' 
#' @return file path to the create bam file `bam_out`
#' @importFrom reticulate import_from_path
#' @export
minimap2_tr_align <- function(mm2_prog_path, fa_file, fq_in, bam_out) {
    callBasilisk(flames_env, function(mm2, fa, fq, bam) {
        python_path <- system.file("python", package="FlamesR")

        align <- reticulate::import_from_path("minimap2_align", python_path)
        align$minimap2_tr_align(mm2, fa, fq, bam)
    }, mm2=mm2_prog_path, fa=fa_file, fq=fq_in, bam=bam_out)

    bam_out # output file
}


