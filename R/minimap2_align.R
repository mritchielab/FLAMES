#' GFF3 to BED12
#'
#' Converts a gff3 file to a bed12, using paftools.js from the minimap2 package
#' 
#' @param minimap2_prog_path Directory containing minimap2, k8 and paftools.js. These are 
#' used to convert gff3 to bed12.
#'
#' @param minimap2_prog_path Absolute path to the directory containing k8, paftools,js and minimap2 
#'      executables from the minimap2 package. Leave as default if these are in the current working directory
#' @param gff3_file The gff3_file to convert
#' @param bed12_file The filename of the bed12 output file.
#' 
#' @importFrom reticulate import_from_path
#' 
#' @export
gff3_to_bed12 <- function(minimap2_prog_path=NULL, gff3_file, bed12_file) {
    python_path <- system.file("python", package="FlamesR")
    cat("\tConverting file", gff3_file, "to bed12:", bed12_file, "\n")
    subprocess_out <- callBasilisk(flames_env, function(mm2_path, gff3, bed12) {
        align <-reticulate::import_from_path("minimap2_align", python_path)
     #   reticulate::source_python(paste(python_path, "minimap2_align.py", sep=.Platform$file.sep))

        if (is.null(mm2_path)) mm2_path = ""
        align$gff3_to_bed12(mm2_path, gff3, bed12)
    }, mm2_path=minimap2_prog_path, gff3=gff3_file, bed12=bed12_file)
    bed12_file # output file
}

#' Minimap2 Align to Genome
#'
#' DESC
#' 
#' @param minimap2_prog_path Absolute path to the directory containing k8, paftools,js and minimap2
#' @param fa_file d
#' @param fq_in d
#' @param bam_out Output BAM file
#' @param no_flank d
#' @param bed12_junc  d
#'
#' @importFrom reticulate import_from_path
#'
#' @export
minimap2_align <- function(minimap2_prog_path=NULL, fa_file, fq_in, bam_out, no_flank=FALSE, bed12_junc=NULL) {
    cat("\tAligning using minimap2. Files:", fa_file, ",", fq_in, "\n")
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
#' Sorts and then indexs alignments for fast random access using samtools.
#' 
#' @param bam_in Input BAM File; MORE
#'
#' @param name Output BAM file; MORE
#' @importFrom reticulate import_from_path
#' @export
samtools_sort_index <- function(bam_in, bam_out) {
    cat("\tWHAT?\n")
    callBasilisk(flames_env, function(bin, bout) {
        python_path <- system.file("python", package="FlamesR")
        sam <- reticulate::import_from_path("minimap2_align", python_path)

        sam$samtools_sort_index(bin, bout)
    }, bin=bam_in, bout=bam_out)
    bam_out
}

