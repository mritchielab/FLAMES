#' Annotate filter (GFF3)
#'
#' @description Combine FLAMES ouput with reference and filter out transcript by
#' ealignment result.
#'
#' @param isoform_gff NEEDED
#' @param ref_gff NEEDED
#' @param isoform_out output isoform file path
#' @param anno_out output annotation file path
#' @param tr_cnt NEEDED
#' @param min_sup_reads NEEDED
#'
#' @importFrom reticulate import_from_path
#' @export
annotate_filter_gff <- function(isoform_gff, ref_gff, isoform_out, anno_out, tr_cnt, min_sup_reads) {
    callBasilisk(flames_env, function(isoform_gff, ref_gff, isoform_out, anno_out, tr_cnt, min_sup_reads) {
        python_path <- system.file("python", package="FlamesR")

        filter <- reticulate::import_from_path("filter_gff", python_path)

        filter$annotate_filter_gff(isoform_gff, ref_gff, isoform_out, anno_out, tr_cnt, min_sup_reads)
    }, isoform_gff=isoform_gff, ref_gff=ref_gff, isoform_out=isoform_out, anno_out=anno_out, tr_cnt=tr_cnt, min_sup_reads=min_sup_reads)

    invisible()
}

