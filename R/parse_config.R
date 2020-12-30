#' Parse Json Configuration File
#'
#' @description Convert a json configuration file into a named R list, grouped into sub lists according to their 
#'      usage in the Flames pipeline.
#' 
#' @param json_file The file name to convert into an R list.
#'
#' @return A named R list of the parameters in \code{json_file}. Subsections are: \code{pipeline_parameters},
#'      \code{global_parameters}, \code{isoform_parameters}, \code{alignment_parameters}, \code{realign_parameters} and
#'      \code{transcript_counting}.
#' @importFrom basilisk basiliskStart basiliskStop basiliskRun
#' @importFrom reticulate import_from_path
#' @export
parse_json_config <- function(json_file) {
    config <- callBasilisk(flames_env, function(json) {
        python_path <- system.file("python", package="FlamesR")

        conf <- reticulate::import_from_path("parse_config", python_path)

        conf$parse_json_config(json)
    }, json=json_file)

    config
}

#' Print Configuration File
#'
#' @details Print the configuration file, represented as a named list used for the Flames pipeline.
#' 
#' @param config List; the configuration list to print.
#'
#' @importFrom reticulate import_from_path
#' @export
print_config <- function(config) {
    callBasilisk(flames_env, function(config) {
        python_path <- system.file("python", package="FlamesR")

        conf <- reticulate::import_from_path("parse_config", python_path)
        conf$print_config(config)
    }, config=config)
    invisible()
}

#' Write Configuration Dictionary to File
#'
#' @details Print the configuration file, represented as a named list used for the Flames pipeline.
#' 
#' @param config List; the configuration list to print.
#' @param config_file the file to output \code{config} to. Should be .json extension
#' @importFrom reticulate import_from_path
write_config <- function(config, config_file) {
    # write the config file to given file path
    callBasilisk(flames_env, function(config, config_file) {
        python_path <- system.file("python", package="FlamesR")

        conf <- reticulate::import_from_path("parse_config", python_path)
        conf$write_config(config, config_file)
    }, config=config, config_file=config_file)
    invisible()
}

#' Write Configuration Dictionary to File
#'
#' @details Print the configuration file, represented as a named list used for the Flames pipeline.
#' 
#' @param do_genome_align
#' @param do_isoform_id
#' @param do_read_realign
#' @param do_transcript_quanti
#' @param gen_raw_isoform
#' @param has_UMI
#' @param MAX_DIST
#' @param MAX_TS_DIST
#' @param MAX_SPLICE_MATCH_DIST
#' @param min_fl_exon_len
#' @param Max_site_per_splice
#' @param Min_sup_cnt
#' @param Min_cnt_pct
#' @param Min_sup_pct
#' @param strand_specific
#' @param remove_incomp_reads
#' @param use_junctions
#' @param no_flank
#' @param use_annotation
#' @param min_tr_coverage
#' @param min_read_coverage
#' @importFrom reticulate import_from_path
create_config <- function(do_genome_align, do_isoform_id,
                            do_read_realign, do_transcript_quanti,
                            gen_raw_isoform, has_UMI,
                            MAX_DIST, MAX_TS_DIST, MAX_SPLICE_MATCH_DIST,
                            min_fl_exon_len, Max_site_per_splice, Min_sup_cnt,
                            Min_cnt_pct, Min_sup_pct, strand_specific, remove_incomp_reads,
                            use_junctions, no_flank,
                            use_annotation, min_tr_coverage, min_read_coverage) {
    # setup config file if none is given, or read in json config
    config =
        list(
            pipeline_parameters=
                list(do_genome_alignment=do_genome_align, do_isoform_identification=do_isoform_id,
                    do_read_realignment=do_read_realign, do_transcript_quantification=do_transcript_quanti),
            global_parameters=
                list(generate_raw_isoform=gen_raw_isoform, has_UMI=has_UMI),
            isoform_parameters=
                list(MAX_DIST=MAX_DIST, MAX_TS_DIST=MAX_TS_DIST, MAX_SPLICE_MATCH_DIST=MAX_SPLICE_MATCH_DIST,
                            min_fl_exon_len=min_fl_exon_len, Max_site_per_splice=Max_site_per_splice, Min_sup_cnt=Min_sup_cnt,
                            Min_cnt_pct=Min_cnt_pct, Min_sup_pct=Min_sup_pct, strand_specific=strand_specific,
                            remove_incomp_reads=remove_incomp_reads),
            alignment_parameters=
                list(use_junctions=use_junctions, no_flank=no_flank),
            realign_parameters=
                list(use_annotation=use_annotation),
            transcript_counting=
                list(min_tr_coverage=min_tr_coverage, min_read_coverage=min_read_coverage)
            )
    if (MAX_DIST <= 0 || MAX_TS_DIST <= 0 || MAX_SPLICE_MATCH_DIST <= 0 || Max_site_per_splice <= 0 ||
        Min_sup_cnt <= 0 || Min_sup_pct <= 0 || (strand_specific != -1 && strand_specific != 0 && strand_specific != 1) || remove_incomp_reads < 0) {
            stop("MAX_DIST,  MAX_TS_DIST, MAX_SPLICE_MATCH_DIST,  Max_site_per_splce, Min_sup_cnt and Min_sup_pct must be greater than 0. strand_specific must be -1, 0 or 1 and remove_incomp_reads must be >= 0.")
    }

    # write created config file.
    config_file_path <- paste(outdir, paste0("config_file_", Sys.getpid(), ".json"), sep="/")
    cat("Writing configuration parameters to: ", config_file_path, "\n")
    write_config(config, config_file_path)

    config
}