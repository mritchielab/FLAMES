#' @importFrom reticulate import_from_path
blocks_to_junctions <- function(block) {
    junctions <- callBasilisk(flames_env, function(block) {
        python_path <- system.file("python", package = "FLAMES")

        sc <-
            reticulate::import_from_path("sc_longread", python_path)
        junc <- sc$blocks_to_junctions(block)
    }, block = block)

    junctions
}

#' @importFrom reticulate import_from_path
remove_similar_tr <-
    function(gene_to_transcript,
             transcript_to_exon,
             thr = 10) {
        callBasilisk(flames_env, function(gene_tran, tr_exon, thr) {
            python_path <- system.file("python", package = "FLAMES")

            sc <-
                reticulate::import_from_path("sc_longread", python_path)
            sc$remove_similar_tr(gene_tran, tr_exon, thr)
        },
        gene_tran = gene_to_transcript, tr_exon = transcript_to_exon, thr =
            thr
        )

        invisible()
    }

#' @importFrom reticulate import_from_path
get_gene_flat <- function(gene_to_transcript, transcript_to_exon) {
    gene_flat <-
        callBasilisk(flames_env, function(gene_tran, tran_exon) {
            python_path <- system.file("python", package = "FLAMES")

            sc <-
                reticulate::import_from_path("sc_longread", python_path)
            sc$get_gene_flat(gene_tran, tran_exon)
        }, gene_tran = gene_to_transcript, tran_exon = transcript_to_exon)

    gene_flat
}

#' @importFrom reticulate import_from_path
get_gene_blocks <-
    function(gene_dict,
             chr_to_gene,
             gene_to_transcript) {
        gene_blocks <-
            callBasilisk(flames_env, function(g_dict, chr_gene, gene_tran) {
                python_path <- system.file("python", package = "FLAMES")

                sc <-
                    reticulate::import_from_path("sc_longread", python_path)

                sc$get_gene_blocks(g_dict, chr_gene, gene_tran)
            }, g_dict = gene_dict, chr_gene = chr_to_gene, gene_tran = gene_to_transcript)

        gene_blocks
    }

#' @importFrom reticulate import_from_path
group_bam2isoform <-
    function(bam_in,
             out_gff3,
             out_stat,
             summary_csv,
             chr_to_blocks,
             gene_dict,
             transcript_to_junctions,
             transcript_dict,
             fa_f,
             config,
             downsample_ratio,
             raw_gff3 = NULL) {
        callBasilisk(flames_env, function(bin,
                                          o_gff3,
                                          o_stat,
                                          summary,
                                          chr,
                                          gene,
                                          trans_junc,
                                          trans_dict,
                                          fa,
                                          conf,
                                          dr,
                                          raw) {
            python_path <- system.file("python", package = "FLAMES")

            sc <-
                reticulate::import_from_path("sc_longread", python_path)
            sc$group_bam2isoform(
                bin,
                o_gff3,
                o_stat,
                summary,
                chr,
                gene,
                trans_junc,
                trans_dict,
                fa,
                conf,
                dr,
                raw
            )
        },
        bin = bam_in, o_gff3 = out_gff3, o_stat = out_stat, summary = summary_csv, chr =
            chr_to_blocks, gene = gene_dict,
        trans_junc = transcript_to_junctions, trans_dict = transcript_dict, fa =
            fa_f, conf = config, dr = downsample_ratio, raw = raw_gff3
        )

        c(out_gff3, out_stat) # output files
    }