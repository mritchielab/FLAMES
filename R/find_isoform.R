#' @importFrom reticulate import_from_path
#' @importFrom Rsamtools indexFa
find_isoform <- function(gff3,
                         genome_bam,
                         isoform_gff3,
                         tss_tes_stat,
                         genomefa,
                         transcript_fa,
                         downsample_ratio,
                         config,
                         raw) {
    ret <- callBasilisk(flames_nopysam_env, function(gff3,
                                                     genome,
                                                     iso,
                                                     tss,
                                                     fa,
                                                     tran,
                                                     ds,
                                                     conf,
                                                     raw) {
        python_path <- system.file("python", package = "FLAMES")

        find <-
            reticulate::import_from_path("find_isoform", python_path)
        ret <-
            find$find_isoform(gff3, genome, iso, tss, fa, tran, ds, conf, raw)

        ret
    },
    gff3 = gff3, genome = genome_bam, iso = isoform_gff3, tss = tss_tes_stat,
    fa = genomefa, tran = transcript_fa, ds = downsample_ratio, conf =
        config, raw = raw
    )

    # we then need to use Rsamtools to index transcript_fa
    Rsamtools::indexFa(transcript_fa) # index the output fa file

    ret
}
