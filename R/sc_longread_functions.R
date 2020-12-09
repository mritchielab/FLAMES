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

#' Title
#'
#' DESC
#' 
#' @param name desc
#'
#' @param name desc
#' @importFrom reticulate import_from_path
#' @export
get_gene_flat <- function(gene_to_transcript, transcript_to_exon) {
    gene_flat <- callBasilisk(flames_env, function(gene_tran, tran_exon) {
        python_path <- system.file("python", package="FlamesR")

        sc <- reticulate::import_from_path("sc_longread", python_path)
        sc$get_gene_flat(gene_tran, tran_exon)
    }, gene_tran=gene_to_transcript, tran_exon=transcript_to_exon)

    gene_flat
}

#' Title
#'
#' DESC
#' 
#' @param name desc
#'
#' @param name desc
#' @importFrom reticulate import_from_path
#' @export
get_gene_blocks <- function(gene_dict, chr_to_gene, gene_to_transcript) {
    gene_blocks <- callBasilisk(flames_env, function(g_dict, chr_gene, gene_tran) {
        python_path <- system.file("python", package="FlamesR")

        sc <- reticulate::import_from_path("sc_longread", python_path)
        
        sc$get_gene_blocks(g_dict, chr_gene, gene_tran)
    }, g_dict=gene_dict, chr_gene=chr_to_gene, gene_tran=gene_to_transcript)

    gene_blocks
}

#' Group BAM to Isoform
#'
#' DESC
#' 
#' @param name desc
#'
#' @param name desc
#' @importFrom reticulate import_from_path
#' @export
group_bam2isoform <- function(bam_in, out_gff3, out_stat, summary_csv, chr_to_blocks, gene_dict,
                             transcript_to_junctions, transcript_dict, fa_f, config, downsample_ratio, 
                             raw_gff3=NULL) {
    # config is a dictionary containing: random_seed, min_cnt_pct
    #["MAX_DIST"]
    #["MAX_TS_DIST"]
    #["MAX_SPLICE_MATCH_DIST"]
    #["Max_site_per_splice"]
    #["Min_sup_cnt"]
    #["min_fl_exon_len"]
    #["Min_sup_pct"]
    #["strand_specific"]
    #["remove_incomp_reads"]
    callBasilisk(flames_env, function(bin, o_gff3, o_stat, summary, chr, gene, trans_junc, trans_dict, fa, conf, dr, raw) {
        python_path <- system.file("python", package="FlamesR")
        
        sc <- reticulate::import_from_path("sc_longread", python_path)
        sc$group_bam2isoform(bin, o_gff3, o_stat, summary, chr, gene, trans_junc, trans_dict, fa, conf, dr, raw)
    }, bin=bam_in, o_gff3=out_gff3, o_stat=out_stat, summary=summary_csv, chr=chr_to_blocks, gene=gene_dict, 
        trans_junc=transcript_to_junctions, trans_dict=transcript_dict, fa=fa_f, conf=config, dr=downsample_ratio, raw=raw_gff3)

    c(out_gff3, out_stat) # output files
}


