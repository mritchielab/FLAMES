#' Parse Json Configuration File
#'
#' @description Convert a json configuration file into a named R list, grouped into sub lists according to their
#'      usage in the Flames pipeline.
#'
#' @param json_file the file path to the JSON file to convert into an R list. 
#' This can be the default FLAMES configuration file found using \code{get_default_config_file()}
#'
#' @return A named R list of the parameters in \code{json_file}. Subsections are: \code{pipeline_parameters},
#'      \code{global_parameters}, \code{isoform_parameters}, \code{alignment_parameters}, \code{realign_parameters} and
#'      \code{transcript_counting}.
#'
#' @examples
#' config <- get_default_config_file()
#' \dontrun{
#' parse_json_config(config)
#' }
#' @importFrom basilisk basiliskStart basiliskStop basiliskRun
#' @importFrom reticulate import_from_path
#' @export
parse_json_config <- function(json_file) {
    config <- callBasilisk(flames_env, function(json) {
        python_path <- system.file("python", package = "FLAMES")

        conf <-
            reticulate::import_from_path("parse_config", python_path)

        conf$parse_json_config(json)
    }, json = json_file)

    config
}

#' Print Configuration File
#'
#' @details Print the configuration file, represented as a named list used for the Flames pipeline.
#'
#' @param config List; the configuration list to print.
#'
#' @return return NULL
#'
#' @importFrom reticulate import_from_path
print_config <- function(config) {
    callBasilisk(flames_env, function(config) {
        python_path <- system.file("python", package = "FLAMES")

        conf <-
            reticulate::import_from_path("parse_config", python_path)
        conf$print_config(config)
    }, config = config)
    invisible()
}

#' Write Configuration Dictionary to File
#'
#' @details Write the configuration file, represented as a named list used for the Flames pipeline.
#'
#' @param config List; the configuration list to write to file.
#' @param config_file the file to output \code{config} to. Should be .json extension
#'
#' @return returns NULL
#'
#' @importFrom reticulate import_from_path
write_config <- function(config, config_file) {
    # write the config file to given file path
    callBasilisk(flames_env, function(config, config_file) {
        python_path <- system.file("python", package = "FLAMES")

        conf <-
            reticulate::import_from_path("parse_config", python_path)
        conf$write_config(config, config_file)
    }, config = config, config_file = config_file)
    invisible()
}

#' Create Configuration File From Arguments
#'
#' @details Create a list object containing the arguments supplied in a format usable for the FLAMES pipeline.
#' Also writes the object to a JSON file, which is located with the prefix 'config_' in the supplied \code{outdir}.
#'
#' @param outdir the destination directory for the configuratio nfile
#' @param do_genome_align Boolean. Specifies whether to run the genome alignment step. \code{TRUE} is recommended
#' @param do_isoform_id Boolean. Specifies whether to run the isoform identification step. \code{TRUE} is recommended
#' @param do_read_realign Boolean. Specifies whether to run the read realignment step. \code{TRUE} is recommended
#' @param do_transcript_quanti Boolean. Specifies whether to run the transcript quantification step. \code{TRUE} is recommended
#' @param gen_raw_isoform Boolean.
#' @param has_UMI Boolean. Specifies if the data contains UMI.
#' @param MAX_DIST Maximum distance allowed when merging splicing sites in isoform consensus clustering.
#' @param MAX_TS_DIST Maximum distance allowed when merging transcript start/end position in isoform consensus clustering.
#' @param MAX_SPLICE_MATCH_DIST Maximum distance allowed when merging splice site called from the data and the reference annotation.
#' @param min_fl_exon_len Minimum length for the first exon outside the gene body in reference annotation. This is to correct the alignment artifact
#' @param Max_site_per_splice Maximum transcript start/end site combinations allowed per splice chain
#' @param Min_sup_cnt Minimum number of read support an isoform decrease this number will significantly increase the number of isoform detected.
#' @param Min_cnt_pct Minimum percentage of count for an isoform relative to total count for the same gene.
#' @param Min_sup_pct Minimum percentage of count for an splice chain that support a given transcript start/end site combination.
#' @param strand_specific 1, -1 or 0. 1 indicates if reads are in the same
#' strand as mRNA, -1 indicates reads are reverse complemented, 0 indicates
#' reads are not strand specific.
#' @param remove_incomp_reads The strenge of truncated isoform filtering. larger number means more stringent filtering.
#' @param use_junctions whether to use known splice junctions to help correct the alignment results
#' @param no_flank Boolean. for synthetic spike-in data. refer to Minimap2 document for detail
#' @param use_annotation Boolean. whether to use reference to help annotate known isoforms
#' @param min_tr_coverage Minimum percentage of isoform coverage for a read to be aligned to that isoform
#' @param min_read_coverage Minimum percentage of read coverage for a read to be uniquely aligned to that isoform
#' @return returns NULL
#' @examples
#' # create the default configuartion file
#' \dontrun{
#' config <- create_config(
#'     getwd(), TRUE, TRUE,
#'     TRUE, TRUE,
#'     TRUE, FALSE,
#'     10, 100, 10,
#'     40, 3, 10,
#'     0.01, 0.2, 1, 5,
#'     TRUE, TRUE,
#'     TRUE, 0.75, 0.75
#' )
#' }
#' @importFrom reticulate import_from_path
#' @export
create_config <- function(outdir,
                          do_genome_align,
                          do_isoform_id,
                          do_read_realign,
                          do_transcript_quanti,
                          gen_raw_isoform,
                          has_UMI,
                          MAX_DIST,
                          MAX_TS_DIST,
                          MAX_SPLICE_MATCH_DIST,
                          min_fl_exon_len,
                          Max_site_per_splice,
                          Min_sup_cnt,
                          Min_cnt_pct,
                          Min_sup_pct,
                          strand_specific,
                          remove_incomp_reads,
                          use_junctions,
                          no_flank,
                          use_annotation,
                          min_tr_coverage,
                          min_read_coverage) {
    # setup config file if none is given, or read in json config
    config <-
        list(
            pipeline_parameters =
                list(
                    do_genome_alignment = do_genome_align,
                    do_isoform_identification = do_isoform_id,
                    do_read_realignment = do_read_realign,
                    do_transcript_quantification = do_transcript_quanti
                ),
            global_parameters =
                list(generate_raw_isoform = gen_raw_isoform, has_UMI = has_UMI),
            isoform_parameters =
                list(
                    MAX_DIST = MAX_DIST,
                    MAX_TS_DIST = MAX_TS_DIST,
                    MAX_SPLICE_MATCH_DIST = MAX_SPLICE_MATCH_DIST,
                    min_fl_exon_len = min_fl_exon_len,
                    Max_site_per_splice = Max_site_per_splice,
                    Min_sup_cnt = Min_sup_cnt,
                    Min_cnt_pct = Min_cnt_pct,
                    Min_sup_pct = Min_sup_pct,
                    strand_specific = strand_specific,
                    remove_incomp_reads = remove_incomp_reads
                ),
            alignment_parameters =
                list(use_junctions = use_junctions, no_flank = no_flank),
            realign_parameters =
                list(use_annotation = use_annotation),
            transcript_counting =
                list(
                    min_tr_coverage = min_tr_coverage,
                    min_read_coverage = min_read_coverage
                )
        )
    if (MAX_DIST <= 0 ||
        MAX_TS_DIST <= 0 ||
        MAX_SPLICE_MATCH_DIST <= 0 || Max_site_per_splice <= 0 ||
        Min_sup_cnt <= 0 ||
        Min_sup_pct <= 0 ||
        (strand_specific != -1 &&
            strand_specific != 0 &&
            strand_specific != 1) || remove_incomp_reads < 0) {
        stop(
            "MAX_DIST,  MAX_TS_DIST, MAX_SPLICE_MATCH_DIST,  Max_site_per_splce, Min_sup_cnt and Min_sup_pct must be greater than 0. strand_specific must be -1, 0 or 1 and remove_incomp_reads must be >= 0."
        )
    }

    # write created config file.
    config_file_path <-
        paste(outdir, paste0("config_file_", Sys.getpid(), ".json"), sep = "/")
    cat(
        "Writing configuration parameters to: ",
        config_file_path,
        "\n"
    )
    write_config(config, config_file_path)

    invisible()
}
