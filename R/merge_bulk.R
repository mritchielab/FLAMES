#' Merge FASTQ (deprecated)
#'
#' Merges all fastq files in the given folder into a single file. Deprecated as the pipeline can
#' now handle multiple fastq files as input.
#'
#' @param fastq_dir Path to the folder containing fastq files to merge
#' @param out_fastq file path to the fastq file which will be created to store all fastq entries. Overwrites existing files
#'
#' @return file path to the created merged fastq file `out_fastq`
#'
#' @examples
#' # download the fastq files to merge
#' temp_path <- tempfile()
#' bfc <- BiocFileCache::BiocFileCache(temp_path, ask = FALSE)
#' file_url <-
#'     "https://raw.githubusercontent.com/OliverVoogd/FLAMESData/master/data"
#' # download the required fastq files, and move them to new folder
#' fastq1 <- bfc[[names(BiocFileCache::bfcadd(bfc, "Fastq1", paste(file_url, "fastq/sample1.fastq.gz", sep = "/")))]]
#' fastq2 <- bfc[[names(BiocFileCache::bfcadd(bfc, "Fastq2", paste(file_url, "fastq/sample2.fastq.gz", sep = "/")))]]
#' fastq_dir <- paste(temp_path, "fastq_dir", sep = "/") # the downloaded fastq files need to be in a directory to be merged together
#' dir.create(fastq_dir)
#' file.copy(c(fastq1, fastq2), fastq_dir)
#' unlink(c(fastq1, fastq2)) # the original files can be deleted
#'
#' # merge the fastq files
#' out_fastq <- paste0(temp_path, "/outfastq.fastq.gz")
#' merge_bulk_fastq(fastq_dir, out_fastq)
#' @importFrom reticulate import_from_path
#' @export
merge_bulk_fastq <- function(fastq_dir, out_fastq) {
    # callBasilisk(flames_nopysam_env, function(fq_dir, out_fq) {
    #     python_path <- system.file("python", package = "FLAMES")

    #     merge_bulk <-
    #         reticulate::import_from_path("merge_bulk_fq", python_path)

    #     merge_bulk$merge_bulk_fq(fq_dir, out_fq)
    # }, fq_dir = fastq_dir, out_fq = out_fastq)

    # out_fastq
    if (utils::file_test("-f", fastq_dir)) {
        merge_bulk_fastq_cpp(fastq_files, out_fastq)
        return(out_fastq)
    }
    fastq_files <- paste(fastq_dir, list.files(fastq_dir), sep = "/")
    fastq_files <- fastq_files[grepl("\\.(fastq|fq)(\\.gz)?$", fastq_files) & utils::file_test("-f", fastq_files)]
    merge_bulk_fastq_cpp(fastq_files, out_fastq)
    return(out_fastq)
}

#' Merge FASTQ using python. (deprecated)
#'
#' Merges all fastq files in the given folder into a single file. Deprecated as the pipeline can
#' now handle multiple fastq files as input.
#'
#' @param fastq_dir Path to the folder containing fastq files to merge
#' @param out_fastq file path to the fastq file which will be created to store all fastq entries. Overwrites existing files
#'
#' @return file path to the created merged fastq file `out_fastq`
#'
#' @examples
#' # download the fastq files to merge
#' temp_path <- tempfile()
#' bfc <- BiocFileCache::BiocFileCache(temp_path, ask = FALSE)
#' file_url <-
#'     "https://raw.githubusercontent.com/OliverVoogd/FLAMESData/master/data"
#' # download the required fastq files, and move them to new folder
#' fastq1 <- bfc[[names(BiocFileCache::bfcadd(bfc, "Fastq1", paste(file_url, "fastq/sample1.fastq.gz", sep = "/")))]]
#' fastq2 <- bfc[[names(BiocFileCache::bfcadd(bfc, "Fastq2", paste(file_url, "fastq/sample2.fastq.gz", sep = "/")))]]
#' fastq_dir <- paste(temp_path, "fastq_dir", sep = "/") # the downloaded fastq files need to be in a directory to be merged together
#' dir.create(fastq_dir)
#' file.copy(c(fastq1, fastq2), fastq_dir)
#' unlink(c(fastq1, fastq2)) # the original files can be deleted
#'
#' # merge the fastq files
#' out_fastq <- paste0(temp_path, "/outfastq.fastq.gz")
#' merge_bulk_fastq(fastq_dir, out_fastq)
#' @importFrom reticulate import_from_path
#' @export
merge_bulk_fastq_python <- function(fastq_dir, out_fastq) {
    cat("WARNING: This function has depreciated. Use FLAMES::merge_bulk_fastq instead.\n")

    callBasilisk(flames_env, function(fq_dir, out_fq) {
        python_path <- system.file("python", package = "FLAMES")

        merge_bulk <-
            reticulate::import_from_path("merge_bulk_fq", python_path)

        merge_bulk$merge_bulk_fq(fq_dir, out_fq)
    }, fq_dir = fastq_dir, out_fq = out_fastq)

    out_fastq
}
