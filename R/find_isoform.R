#' Title
#'
#' DESC
#' 
#' @param name desc
#'
#' @param name desc
#' @importFrom reticulate import_from_path
#' @export
find_isoform <- function(anno, genome_bam, isoform_gff3, tss_tes_stat, genomefa, do_isoform_identification) {
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
        group_bam2isoform(genome_bam, isoform_gff3, tss_tes_stat, "", chr_to_blocks, 
                gene_dict, transcript_to_junctions, parse_res$transcript_dict, genomefa,
                config=isoform_parameters, downsample_ratio=downsample_ratio, 
                raw_gff3=if (generate_raw_isoform) raw_splice_isoform else NULL)
    } else {
        ## skip finding isoform.
    }
}