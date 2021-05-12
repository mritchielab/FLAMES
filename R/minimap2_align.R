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
gff3_to_bed12 <- function(minimap2_prog_path=NULL, gff3_file, bed12_file) {
    python_path <- system.file("python", package="FLAMES")
    callBasilisk(flames_nopysam_env, function(mm2_path, gff3, bed12) {
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
#' Uses options "-ax splice -t 12 -k14 --secondary=no \code{fa_file} \code{fq_in}"
#' 
#' @param minimap2_prog_path Absolute path to the directory containing minimap2
#' @param fa_file Fasta file used as a reference database for alignment
#' @param fq_in Fastq file used as a query sequence file
#' @param sam_out Output SAM file
#' @param no_flank Boolean; used if studying SIRV, to let minimap2 ignore additional bases
#' @param bed12_junc Gene annotations in BED12 format. If specified, minmap2 prefers splicing in annotations.
#'
#' @return file path to the given output BAM file, \code{bam_out}
#' @importFrom reticulate import_from_path
minimap2_align <- function(minimap2_prog_path=NULL, fa_file, fq_in, sam_out, no_flank=FALSE, bed12_junc=NULL) {
    callBasilisk(flames_nopysam_env, function (mm2_path, fa, fq, sam, flank, bed12_junc) {
        python_path <- system.file("python", package="FLAMES")
        mm2 <- reticulate::import_from_path("minimap2_align", python_path)

        if (is.null(minimap2_prog_path)) minimap2_prog_path = ""
        mm2$minimap2_align(mm2_path, fa, fq, sam, flank, bed12_junc)
    }, mm2_path=minimap2_prog_path, fa=fa_file, fq=fq_in, sam=sam_out, flank=no_flank, bed12_junc=bed12_junc)
    sam_out # output file
}

#' Samtools Sort and Index
#' 
#' @description 
#' Sort and index the given BAM file, using Rsamtools.
#' 
#' @param bam_in the bam file to sort and index
#' @param bam_out the output indexed bam file.
#' @return the same value as \code{bam_out}
#' @importFrom Rsamtools sortBam indexBam
#' 
#' @examples 
#' temp_path <- tempfile()
#' bfc <- BiocFileCache::BiocFileCache(temp_path, ask=FALSE)
#' file_url <- 
#'    "https://raw.githubusercontent.com/OliverVoogd/FLAMESData/master/data"
#' genome_bam <- paste0(temp_path, "/align2genome.bam")
#' file.rename(bfc[[names(BiocFileCache::bfcadd(bfc, "Genome BAM", paste(file_url, "align2genome.bam", sep="/")))]], genome_bam)
#' 
#' samtools_sort_index(genome_bam, tempfile(fileext=".bam"))
#' @export
samtools_sort_index <- function(bam_in, bam_out) {
    Rsamtools::sortBam(bam_in, gsub("\\.bam", "", bam_out))
    Rsamtools::indexBam(bam_out)

    bam_out
}

#' Samtools as BAM
#' 
#' @description 
#' Produces a compressed binary BAM file from a text based SAM file, using Rsamtools.
#' 
#' @param sam_in the input SAM file
#' @param bam_out the output BAM file
#' 
#' @return the same value as \code{bam_out}
#' 
#' @examples 
#' temp_path <- tempfile()
#' bfc <- BiocFileCache::BiocFileCache(temp_path, ask=FALSE)
#' file_url <- 
#'    "https://raw.githubusercontent.com/OliverVoogd/FLAMESData/master/data"
#' tmp_sam <- paste0(temp_path, "/tmp_sam.sam")
#' file.rename(bfc[[names(BiocFileCache::bfcadd(bfc, "Temp SAM", paste(file_url, "temp_sam.sam", sep="/")))]], tmp_sam)
#' 
#' samtools_as_bam(tmp_sam, tempfile(fileext="bam"))
#' @importFrom Rsamtools asBam
#' @export
samtools_as_bam <- function(sam_in, bam_out) {
    Rsamtools::asBam(sam_in, gsub("\\.bam", "", bam_out))
    bam_out
}

#' @importFrom reticulate import_from_path
minimap2_tr_align <- function(mm2_prog_path, fa_file, fq_in, sam_out) {
    callBasilisk(flames_nopysam_env, function(mm2, fa, fq, sam) {
        python_path <- system.file("python", package="FLAMES")

        align <- reticulate::import_from_path("minimap2_align", python_path)
        align$minimap2_tr_align(mm2, fa, fq, sam)
    }, mm2=mm2_prog_path, fa=fa_file, fq=fq_in, sam=sam_out)

    sam_out # output file
}


