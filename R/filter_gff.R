#' @importFrom reticulate import_from_path
#' @importFrom basilisk basiliskRun
annotate_filter_gff <-
    function(isoform_gff,
             ref_gff,
             isoform_out,
             anno_out,
             tr_cnt,
             min_sup_reads) {
        basiliskRun(env = flames_env, fun = function(isoform_gff,
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


#' @importFrom reticulate import_from_path
#' @importFrom basilisk basiliskRun
annotate_full_splice_match_all_sample <-
    function(anno_out, isoform_gff, ref_gff) {
        basiliskRun(env = flames_env, fun = function(anno_out, isoform_gff, ref_gff) {
            python_path <- system.file("python", package = "FLAMES")

            filter <-
                reticulate::import_from_path("filter_gff", python_path)

            filter$annotate_full_splice_match_all_sample(anno_out, isoform_gff, ref_gff)
        },
        anno_out = anno_out, isoform_gff = isoform_gff, ref_gff = ref_gff
        )

        invisible()
    }
