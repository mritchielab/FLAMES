#' Demultiplex reads using Sockeye outputs 
#' @description Demultiplex reads using the \code{cell_umi_gene.tsv} file from Sockeye.
#' @param fastq_dir The folder containing FASTQ files from Sockeye's output under \code{ingest/chunked_fastqs}.
#' @param sockeye_tsv The \code{cell_umi_gene.tsv} file from Sockeye.
#' @param out_fq The output FASTQ file.
#' @return NULL
#' @export
#' @importFrom reticulate import_from_path
demultiplex_sockeye <- function(fastq_dir, sockeye_tsv, out_fq) {
    callBasilisk(flames_env, function(fastq_dir_, sockeye_tsv_, out_fq_) {
        python_path <- system.file("python", package = "FLAMES")
        demult <- reticulate::import_from_path("demultiplex_sockeye", python_path)
        ret <- demult$demultiplex_sockeye(fastq_dir_, sockeye_tsv_, out_fq_)
        ret
    }, fastq_dir_ = fastq_dir, sockeye_tsv_ = sockeye_tsv, out_fq_ = out_fq)
}
