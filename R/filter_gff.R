#' @importFrom reticulate import_from_path
annotate_filter_gff <-
    function(isoform_gff,
             ref_gff,
             isoform_out,
             anno_out,
             tr_cnt,
             min_sup_reads) {
        callBasilisk(flames_nopysam_env, function(isoform_gff,
                                                  ref_gff,
                                                  isoform_out,
                                                  anno_out,
                                                  tr_cnt,
                                                  min_sup_reads) {
            python_path <- system.file("python", package = "FLAMES")

            filter <-
                reticulate::import_from_path("filter_gff", python_path)

            filter$annotate_filter_gff(
                isoform_gff,
                ref_gff,
                isoform_out,
                anno_out,
                tr_cnt,
                min_sup_reads
            )
        },
        isoform_gff = isoform_gff, ref_gff = ref_gff, isoform_out = isoform_out, anno_out =
            anno_out, tr_cnt = tr_cnt, min_sup_reads = min_sup_reads
        )

        invisible()
    }
