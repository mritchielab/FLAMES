#' @importFrom reticulate import_from_path
#' @importFrom Rsamtools indexFa
get_transcript_seq <- function(fa_file, fa_out_f, chr_to_gene, transcript_dict,
                                gene_to_transcript, transcript_to_exon, ref_dict=NULL) {
    callBasilisk(flames_env, function(fa_file, fa_out_f, chr_to_gene, transcript_dict,
                       gene_to_transcript, transcript_to_exon, ref_dict) {
        python_path <- system.file("python", package="FLAMES")

        g_fa <- reticulate::import_from_path("gff3_to_fa", python_path)
        g_fa$get_transcript_seq(fa_file, fa_out_f, chr_to_gene, transcript_dict,
                       gene_to_transcript, transcript_to_exon, ref_dict)
    }, fa_file=fa_file, fa_out_f=fa_out_f, chr_to_gene=chr_to_gene, transcript_dict=transcript_dict,
                       gene_to_transcript=gene_to_transcript, transcript_to_exon=transcript_to_exon, ref_dict=ref_dict)
    cat("TEMP: finished get_transcript_seq, now indexing")
    #Rsamtools::indexFa(file=fa_out_f) # index the output fa file
    invisible()
}

