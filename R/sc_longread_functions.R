#' Title
#'
#' DESC
#' 
#' @param name desc
#'
#' @param name desc
#' @importFrom reticulate import_from_path
#' @export
blocks_to_junctions <- function(block) {
    junctions <- callBasilisk(flames_env, function(block) {
        python_path <- system.file("python", package="FlamesR")

        sc <- reticulate::import_from_path("sc_longread", python_path)
        junc <- sc$blocks_to_junctions(block)
    }, block=block)

    junctions
}

#' Remove Similar Transcripts
#'
#' DESC
#' 
#' @param name desc
#'
#' @param name desc
#' @importFrom reticulate import_from_path
#' @export
remove_similar_tr <- function(transcript_dict, gene_to_transcript, transcript_to_exon, thr=10) {
    callBasilisk(flames_env, function(tr_dict, gene_tran, tr_exon, thr) {
        python_path <- system.file("python", package="FlamesR")

        sc <- reticulate::import_from_path("sc_longread", python_path)
        sc$remove_similar_tr(tr_dict, gene_tran, tr_exon, thr)
    }, tr_dict=transcript_dict, gene_tran=gene_to_transcript, tr_exon=transcript_to_exon, thr=thr)

    invisible()
}

