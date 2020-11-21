#' Title
#'
#' DESC
#' 
#' @param name desc
#'
#' @param name desc
#' @importFrom reticulate import_from_path
#' @export
find_isoform <- function(anno, do_isoform_identification) {
    parse_res <- parse_gff_tree(anno)
    transcript_to_junctions = list()
    for (tr in names(parse_res$transcript_to_exon)) {
        transcript_to_junctions[[tr]] = blocks_to_junctions(parse_res$transcript_to_exon[[tr]])
    }

    remove_similar_tr(parse_res$transcript_dict, parse_res$gene_to_transcript, parse_res$transcript_to_exon)
    gene_dict <- get_gene_flat(parse_res$gene_to_transcript, parse_res$transcript_to_exon)
    chr_to_blocks <- get_gene_blocks(gene_dict, parse_res$chr_to_gene, parse_res$gene_to_transcript)

    if (do_isoform_identification) {
        cat("\nFind isoform")
    }
}