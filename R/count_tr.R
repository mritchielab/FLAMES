#' Parse Realigned Bam
#' @returns realigned bam
#' @examples 
#' 
#' realign_bam <- system.file("extdata/align2genome.bam", package="FLAMES")
#' transcript_fa <- system.file("extdata/transcript_assembly.fa.fai", package="FLAMES")
#' \dontrun{
#' parse_realign <- parse_realigned_bam(realign_bam, transcript_fa, 10, 0.75, 0.75)
#'}
#' @importFrom reticulate import_from_path dict
parse_realigned_bam <- function(bam_in, fa_idx_f, min_sup_reads, min_tr_coverage, min_read_coverage, ...) {

    kwargs = reticulate::dict(...)
    ret <- callBasilisk(flames_env, function(bam_in, fa_idx_f, min_sup_reads, min_tr_coverage, min_read_coverage, kwargs) {
        python_path <- system.file("python", package="FLAMES")

        count <- reticulate::import_from_path("count_tr", python_path)
        ret <- count$parse_realigned_bam(bam_in, fa_idx_f, min_sup_reads, min_tr_coverage, min_read_coverage, kwargs)

        names(ret) <- c("bc_tr_count_dict", "bc_tr_badcov_count_dict", "tr_kept")
        ret
    }, bam_in=bam_in, fa_idx_f=fa_idx_f, min_sup_reads=min_sup_reads, min_tr_coverage=min_tr_coverage, min_read_coverage=min_read_coverage, kwargs=kwargs)

    ret
}

#' Write Transcript to CSV file
#' @return returns NULL
#' @examples 
#' isoform_gff3 <- parse_gff_tree(system.file("extdata/isoform_annotated.gff3", package="FLAMES"))
#' gff3_parse <- parse_gff_tree(system.file("extdata/SIRV_anno.gtf", package="FLAMES"))
#' realign_bam <- system.file("extdata/align2genome.bam", package="FLAMES")
#' \dontrun{
#' parse_realign <- parse_realigned_bam(realign_bam, system.file("extdata/transcript_assembly.fa.fai", package="FLAMES"), 10, 0.75, 0.75)
#' tr_cnt <- wrt_tr_to_csv(parse_realign$bc_tr_count_dict, isoform_gff3$transcript_dict, tempfile(fileext=".csv.gz"), 
#'                         gff3_parse$transcript_dict, FALSE)
#'                         }
#' @importFrom reticulate import_from_path
wrt_tr_to_csv <- function(bc_tr_count_dict, transcript_dict, csv_f, transcript_dict_ref=NULL, has_UMI=TRUE) {
    callBasilisk(flames_env, function(bc_tr_count_dict, transcript_dict, csv_f, transcript_dict_ref, has_UMI) {
        python_path <- system.file("python", package="FLAMES")

        count <- reticulate::import_from_path("count_tr", python_path)

        count$wrt_tr_to_csv(bc_tr_count_dict, transcript_dict, csv_f, transcript_dict_ref, has_UMI)
    }, bc_tr_count_dict=bc_tr_count_dict, transcript_dict=transcript_dict, csv_f=csv_f, transcript_dict_ref=transcript_dict_ref, has_UMI=has_UMI)
}



