#' Blocks To Junctions
#'
#' @description
#' Converts \code{block} to a named list containing the keys "left", "right"
#' and "junctions".
#' @return The converted junctions named list
#' 
#' @examples 
#' # parse the gff file before creating junctions from the gene blocks
#' gff3_parse <- parse_gff_tree(system.file("extdata/SIRV_anno.gtf", package="FLAMES"))
#' transcript_junction = blocks_to_junctions(gff3_parse$transcript_to_exon[[1]])
#' @importFrom reticulate import_from_path
blocks_to_junctions <- function(block) {
    junctions <- callBasilisk(flames_env, function(block) {
        python_path <- system.file("python", package="FLAMES")

        sc <- reticulate::import_from_path("sc_longread", python_path)
        junc <- sc$blocks_to_junctions(block)
    }, block=block)

    junctions
}

#' Remove Similar Transcripts
#'
#' @description
#' Remove any duplicate exons, or similar exons with a similarity greater than
#' \code{thr}, from each element in \code{gene_to_transcript}. Modifies
#' \code{gene_to_transcript} in order to remove duplicates.
#' @return returns NULL
#' 
#' @examples 
#' gff3_parse <- parse_gff_tree(system.file("extdata/SIRV_anno.gtf", package="FLAMES"))
#' remove_similar_tr(gff3_parse$gene_to_transcript, gff3_parse$transcript_to_exon, thr=5)
#' @importFrom reticulate import_from_path
remove_similar_tr <- function(gene_to_transcript, transcript_to_exon, thr=10) {
    callBasilisk(flames_env, function(gene_tran, tr_exon, thr) {
        python_path <- system.file("python", package="FLAMES")

        sc <- reticulate::import_from_path("sc_longread", python_path)
        sc$remove_similar_tr(gene_tran, tr_exon, thr)
    }, gene_tran=gene_to_transcript, tr_exon=transcript_to_exon, thr=thr)

    invisible()
}

#' Get Gene Flat
#' @return gene_flat
#' @examples 
#' gff3_parse <- parse_gff_tree(system.file("extdata/SIRV_anno.gtf", package="FLAMES"))
#' gene_dict <- get_gene_flat(gff3_parse$gene_to_transcript, gff3_parse$transcript_to_exon)
#' @importFrom reticulate import_from_path
get_gene_flat <- function(gene_to_transcript, transcript_to_exon) {
    gene_flat <- callBasilisk(flames_env, function(gene_tran, tran_exon) {
        python_path <- system.file("python", package="FLAMES")

        sc <- reticulate::import_from_path("sc_longread", python_path)
        sc$get_gene_flat(gene_tran, tran_exon)
    }, gene_tran=gene_to_transcript, tran_exon=transcript_to_exon)

    gene_flat
}

#' Get Gene Blocks
#' @examples 
#' gff3_parse <- parse_gff_tree(system.file("extdata/SIRV_anno.gtf", package="FLAMES"))
#' gene_dict <- get_gene_flat(gff3_parse$gene_to_transcript, gff3_parse$transcript_to_exon)
#' chr_to_blocks <- get_gene_blocks(gene_dict, gff3_parse$chr_to_gene, gff3_parse$gene_to_transcript)
#' @importFrom reticulate import_from_path
get_gene_blocks <- function(gene_dict, chr_to_gene, gene_to_transcript) {
    gene_blocks <- callBasilisk(flames_env, function(g_dict, chr_gene, gene_tran) {
        python_path <- system.file("python", package="FLAMES")

        sc <- reticulate::import_from_path("sc_longread", python_path)

        sc$get_gene_blocks(g_dict, chr_gene, gene_tran)
    }, g_dict=gene_dict, chr_gene=chr_to_gene, gene_tran=gene_to_transcript)

    gene_blocks
}

#' Group BAM to Isoform
#' @details \code{config} contains
#' \itemize{
#'     \item MAX_DIST
#'     \item MAX_TS_DIST
#'     \item MAX_SPLICE_MATCH_DIST
#'     \item Max_site_per_splice
#'     \item Min_sup_cnt
#'     \item min_fl_exon_len
#'     \item Min_sup_pct
#'     \item strand_specific
#'     \item remove_incomp_reads
#'     \item random_seed OPTIONAL;
#' }
#' @return File paths of the output files \code{out_gff3} and \code{out_stat}
#' 
#' @examples 
#' anno <- system.file("extdata/SIRV_anno.gtf", package="FLAMES")
#' genome_bam <- system.file("extdata/align2genome.bam", package="FLAMES")
#' genomefa <- system.file("extdata/SIRV_genomefa.fasta", package="FLAMES")
#' config <- parse_json_config(system.file("extdata/SIRV_config_default.json", package="FLAMES"))
#' gff3_parse <- parse_gff_tree(anno)
#' gene_dict <- get_gene_flat(gff3_parse$gene_to_transcript, gff3_parse$transcript_to_exon)
#' chr_to_blocks <- get_gene_blocks(gene_dict, gff3_parse$chr_to_gene, gff3_parse$gene_to_transcript)
#' 
#' juncs = list()
#' for (tr in names(gff3_parse$transcript_to_exon)) juncs[[tr]] = blocks_to_junctions(gff3_parse$transcript_to_exon[[tr]])
#' group_bam2isoform(genome_bam, out_gff3=tempfile(fileext=".gff3"), out_stat=tempfile(fileext=".bedgraph"),
#'                  "", chr_to_blocks, gene_dict, juncs, gff3_parse$transcript_dict, 
#'                  genomefa, config$isoform_parameters, 1, NULL)
#' @importFrom reticulate import_from_path
group_bam2isoform <- function(bam_in, out_gff3, out_stat, summary_csv, chr_to_blocks, gene_dict,
                             transcript_to_junctions, transcript_dict, fa_f, config, downsample_ratio,
                             raw_gff3=NULL) {
    callBasilisk(flames_env, function(bin, o_gff3, o_stat, summary, chr, gene, trans_junc, trans_dict, fa, conf, dr, raw) {
        python_path <- system.file("python", package="FLAMES")

        sc <- reticulate::import_from_path("sc_longread", python_path)
        sc$group_bam2isoform(bin, o_gff3, o_stat, summary, chr, gene, trans_junc, trans_dict, fa, conf, dr, raw)
    }, bin=bam_in, o_gff3=out_gff3, o_stat=out_stat, summary=summary_csv, chr=chr_to_blocks, gene=gene_dict,
        trans_junc=transcript_to_junctions, trans_dict=transcript_dict, fa=fa_f, conf=config, dr=downsample_ratio, raw=raw_gff3)

    c(out_gff3, out_stat) # output files
}


