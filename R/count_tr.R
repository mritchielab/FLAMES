#' Title
#'
#' DESC
#' 
#' @param name desc
#'
#' @param name desc
#' @importFrom reticulate import_from_path dict
#' @export
parse_realigned_bam <- function(bam_in, fa_idx_f, min_sup_reads, min_tr_coverage, min_read_coverage, ...) {
    
    kwargs = reticulate::dict(...)
    ret <- callBasilisk(flames_env, function(bam_in, fa_idx_f, min_sup_reads, min_tr_coverage, min_read_coverage, kwargs) {
        python_path <- system.file("python", package="FlamesR")

        count <- reticulate::import_from_path("count_tr", python_path)
        ret <- count$parse_realigned_bam(bam_in, fa_idx_f, min_sup_reads, min_tr_coverage, min_read_coverage, kwargs)
    
        names(ret) <- c("bc_tr_count_dict", "bc_tr_badcov_count_dict", "tr_kept")
        ret
    }, bam_in=bam_in, fa_idx_f=fa_idx_f, min_sup_reads=min_sup_reads, min_tr_coverage=min_tr_coverage, min_read_coverage=min_read_coverage, kwargs=kwargs)
    
    ret
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
wrt_tr_to_csv <- function(bc_tr_count_dict, transcript_dict, csv_f, transcript_dict_ref=NULL, has_UMI=TRUE) {
    callBasilisk(flames_env, function(bc_tr_count_dict, transcript_dict, csv_f, transcript_dict_ref, has_UMI) {
        python_path <- system.file("python", package="FlamesR")

        count <- reticulate::import_from_path("count_tr", python_path)

        count$wrt_tr_to_csv(bc_tr_count_dict, transcript_dict, csv_f, transcript_dict_ref, has_UMI)
    }, bc_tr_count_dict=bc_tr_count_dict, transcript_dict=transcript_dict, csv_f=csv_f, transcript_dict_ref=transcript_dict_ref, has_UMI=has_UMI)
}



