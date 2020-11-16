#' Title
#'
#' Merges all fastq files in the given folder into a single file
#'
#'  THIS NEED PROPER TYPES
#' @param fastq_dir The folder containing fastq files to merge
#' @param anno_csv a
#' @param out_fastq A file which will be created to store all fastq entries
#' @importFrom reticulate import_from_path
#' @export
merge_bulk_fastq <- function(fastq_dir, anno_csv, out_fastq) {
    callBasilisk(flames_env, function(fq_dir, a_csv, out_fq) {
        python_path <- system.file("python", package="FlamesR")
        
        merge_bulk <-
            reticulate::import_from_path("merge_bulk_fq", python_path)
    
        merge_bulk$merge_bulk_fq(fq_dir, a_csv, out_fq)
        
    }, fq_dir=fastq_dir, a_csv=anno_csv, out_fq=out_fastq)

    out_fastq
}
