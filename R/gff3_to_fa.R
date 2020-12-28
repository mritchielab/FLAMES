#' Get Transcript Sequences
#'
#' @description NEEDED
#'
#' @param fa_file NEEDER
#' @param fa_out_f NEEDED
#' @param chr_to_gene NEEDED
#' @param transcript_dict NEEDED
#' @param gene_to_transcript NEEDED
#' @param transcript_to_exon NEEDED
#' @param ref_dict NEEDED
#' @importFrom reticulate import_from_path
#' @export
get_transcript_seq <- function(fa_file, fa_out_f, chr_to_gene, transcript_dict,
                                gene_to_transcript, transcript_to_exon, ref_dict=NULL) {
    callBasilisk(flames_env, function(fa_file, fa_out_f, chr_to_gene, transcript_dict,
                       gene_to_transcript, transcript_to_exon, ref_dict) {
        python_path <- system.file("python", package="FlamesR")

        g_fa <- reticulate::import_from_path("gff3_to_fa", python_path)
        g_fa$get_transcript_seq(fa_file, fa_out_f, chr_to_gene, transcript_dict,
                       gene_to_transcript, transcript_to_exon, ref_dict)
    }, fa_file=fa_file, fa_out_f=fa_out_f, chr_to_gene=chr_to_gene, transcript_dict=transcript_dict,
                       gene_to_transcript=gene_to_transcript, transcript_to_exon=transcript_to_exon, ref_dict=ref_dict)

    invisible()
}

