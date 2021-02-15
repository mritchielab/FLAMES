#' Merge FASTQ
#'
#' Merges all fastq files in the given folder into a single file. Used to create a pseudobulk fastq file,
#' with added fake barcodes to differentiate between source files.
#'
#' @param fastq_dir The folder containing fastq files to merge
#' @param anno_csv a path for the output csv, containing the fake barcodes created
#' @param out_fastq A file which will be created to store all fastq entries
#' 
#' @return file path to the created merged fastq file `out_fastq`
#' 
#' @examples
#' out_fastq <- merge_bulk_fastq(system.file("extdata/fastq", package="FlamesR"), anno_csv=tempfile(), out_fastq=tempfile())
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
