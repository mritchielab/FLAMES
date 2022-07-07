#' @importFrom reticulate import_from_path dict
parse_realigned_bam <-
    function(bam_in,
             fa_idx_f,
             min_sup_reads,
             min_tr_coverage,
             min_read_coverage,
             bc_file) {
        ret <-
            callBasilisk(flames_env, function(bam_in,
                                              fa_idx_f,
                                              min_sup_reads,
                                              min_tr_coverage,
                                              min_read_coverage,
                                              bc_file) {
                python_path <- system.file("python", package = "FLAMES")

                count <-
                    reticulate::import_from_path("count_tr", python_path)
                ret <-
                    count$parse_realigned_bam(
                        bam_in,
                        fa_idx_f,
                        min_sup_reads,
                        min_tr_coverage,
                        min_read_coverage,
                        bc_file
                    )

                names(ret) <-
                    c(
                        "bc_tr_count_dict",
                        "bc_tr_badcov_count_dict",
                        "tr_kept"
                    )
                ret
            },
            bam_in = bam_in, fa_idx_f = fa_idx_f, min_sup_reads = min_sup_reads, min_tr_coverage =
                min_tr_coverage, min_read_coverage = min_read_coverage, bc_file = bc_file
            )

        ret
    }

#' @importFrom reticulate import_from_path
wrt_tr_to_csv <-
    function(bc_tr_count_dict,
             transcript_dict,
             csv_f,
             transcript_dict_ref = NULL,
             has_UMI = TRUE) {
        callBasilisk(flames_env, function(bc_tr_count_dict,
                                          transcript_dict,
                                          csv_f,
                                          transcript_dict_ref,
                                          has_UMI) {
            python_path <- system.file("python", package = "FLAMES")

            count <-
                reticulate::import_from_path("count_tr", python_path)

            count$wrt_tr_to_csv(
                bc_tr_count_dict,
                transcript_dict,
                csv_f,
                transcript_dict_ref,
                has_UMI
            )
        },
        bc_tr_count_dict = bc_tr_count_dict, transcript_dict = transcript_dict, csv_f =
            csv_f, transcript_dict_ref = transcript_dict_ref, has_UMI = has_UMI
        )
    }

#' @importFrom reticulate import_from_path dict
flames_quantify <- function(annotation, outdir, bulk, config) {
    if (missing("bulk")) {
        if (file_test("-f", file.path(outdir, "pseudo_barcode_annotation.csv"))) {
            bulk <- TRUE
        } else {
            bulk <- FALSE
        }
    }

    callBasilisk(flames_env, function(config_dict, realign_bam, transcript_fa_idx, tr_cnt_csv, isoform_gff3, annotation, isoform_gff3_f, FSM_anno_out, tr_badcov_cnt_csv, bc_file) {
        python_path <- system.file("python", package = "FLAMES")
        count <- reticulate::import_from_path("count_tr", python_path)
        count$quantification(config_dict, realign_bam, transcript_fa_idx, tr_cnt_csv, isoform_gff3, annotation, isoform_gff3_f, FSM_anno_out, tr_badcov_cnt_csv, bc_file)
    },
    config_dict = reticulate::dict(config),
    realign_bam = file.path(outdir, "realign2transcript.bam"),
    transcript_fa_idx = file.path(outdir, "transcript_assembly.fa.fai"),
    tr_cnt_csv = file.path(outdir, "transcript_count.csv.gz"),
    isoform_gff3 = file.path(outdir, ifelse(config$pipeline_parameters$bambu_isoform_identification, "isoform_annotated.gtf", "isoform_annotated.gff3")),
    annotation = annotation,
    isoform_gff3_f = file.path(outdir, "isoform_annotated.filtered.gff3"),
    FSM_anno_out = file.path(outdir, "isoform_FSM_annotation.csv"),
    tr_badcov_cnt_csv = file.path(outdir, "transcript_count.bad_coverage.csv.gz"),
    bc_file = ifelse(bulk, file.path(outdir, "pseudo_barcode_annotation.csv"), FALSE)
    )
}
