#' Create Configuration File From Arguments
#'
#' @details Create a list object containing the arguments supplied in a format usable for the FLAMES pipeline.
#' Also writes the object to a JSON file, which is located with the prefix 'config_' in the supplied \code{outdir}.
#' Default values from \code{extdata/config_sclr_nanopore_3end.json} will be used for unprovided parameters.
#'
#' @param outdir the destination directory for the configuratio nfile
#' @param type use an example config, available values:
#' \itemize{
#'  \item{"sc_5end"}{ - config for 5' end ONT reads}
#'  \item{"SIRV"}{ - config for the SIRV example reads}
#' }
#' @param ... Configuration parameters.
#' \itemize{
#'  \item{do_genome_align}{ - Boolean. Specifies whether to run the genome alignment step. \code{TRUE} is recommended}
#'  \item{do_isoform_id}{ - Boolean. Specifies whether to run the isoform identification step. \code{TRUE} is recommended}
#'  \item{do_read_realign}{ - Boolean. Specifies whether to run the read realignment step. \code{TRUE} is recommended}
#'  \item{do_transcript_quanti}{ - Boolean. Specifies whether to run the transcript quantification step. \code{TRUE} is recommended}
#'  \item{gen_raw_isoform}{ - Boolean.}
#'  \item{has_UMI}{ - Boolean. Specifies if the data contains UMI.}
#'  \item{max_dist}{ - Maximum distance allowed when merging splicing sites in isoform consensus clustering.}
#'  \item{max_ts_dist}{ - Maximum distance allowed when merging transcript start/end position in isoform consensus clustering.}
#'  \item{max_splice_match_dist}{ - Maximum distance allowed when merging splice site called from the data and the reference annotation.}
#'  \item{min_fl_exon_len}{ - Minimum length for the first exon outside the gene body in reference annotation. This is to correct the alignment artifact}
#'  \item{max_site_per_splice}{ - Maximum transcript start/end site combinations allowed per splice chain}
#'  \item{min_sup_cnt}{ - Minimum number of read support an isoform decrease this number will significantly increase the number of isoform detected.}
#'  \item{min_cnt_pct}{ - Minimum percentage of count for an isoform relative to total count for the same gene.}
#'  \item{min_sup_pct}{ - Minimum percentage of count for an splice chain that support a given transcript start/end site combination.}
#'  \item{strand_specific}{ - 0, 1 or -1. 1 indicates if reads are in the same strand as mRNA, -1 indicates reads are reverse complemented, 0 indicates reads are not strand specific.}
#'  \item{remove_incomp_reads}{ - The strenge of truncated isoform filtering. larger number means more stringent filtering.}
#'  \item{use_junctions}{ - whether to use known splice junctions to help correct the alignment results}
#'  \item{no_flank}{ - Boolean. for synthetic spike-in data. refer to Minimap2 document for detail}
#'  \item{use_annotation}{ - Boolean. whether to use reference to help annotate known isoforms}
#'  \item{min_tr_coverage}{ - Minimum percentage of isoform coverage for a read to be aligned to that isoform}
#'  \item{min_read_coverage}{ - Minimum percentage of read coverage for a read to be uniquely aligned to that isoform}
#' }
#'
#' @return file path to the config file created
#' @examples
#' # create the default configuration file
#' outdir <- tempdir()
#' config <- create_config(outdir)
#'
#' @importFrom jsonlite toJSON fromJSON
#' @importFrom utils modifyList
#' @importFrom stats setNames
#' @export
create_config <- function(outdir, type = "sc_3end", ...) {
    if (type == "sc_3end") {
        config <- jsonlite::fromJSON(system.file("extdata/config_sclr_nanopore_3end.json", package = "FLAMES"))
    } else if (type == "SIRV") {
        config <- jsonlite::fromJSON(system.file("extdata/SIRV_config_default.json", package = "FLAMES"))
    } else {
        stop("Unrecognised config type ", type)
    }
    updates <- list(...)

    if (length(updates) > 0) {
        if (any(is.null(names(updates))) || "" %in% names(updates)) {
            stop("Parameters must be named")
        }
        config <- within(config, rm(comment))
        for (i_param in names(updates)) {
            i_part <- names(config)[as.logical(lapply(lapply(config, names), function(part) {
                i_param %in% part
            }))]
            config <- modifyList(config, setNames(list(updates[i_param]), i_part))
        }
    }

    # write created config file.
    config_file_path <- file.path(outdir, paste0("config_file_", Sys.getpid(), ".json"))
    cat(
        "Writing configuration parameters to: ",
        config_file_path,
        "\n"
    )
    write(jsonlite::toJSON(config, pretty = TRUE), config_file_path)

    return(config_file_path)
}

#' @importFrom Matrix tail
#' @importFrom stringr str_split
#' @importFrom jsonlite fromJSON
check_arguments <-
    function(annotation,
             fastq,
             genome_bam,
             outdir,
             genome_fa,
             minimap2_dir,
             config_file) {
        if (!dir.exists(outdir)) {
            cat("Output directory does not exists: one is being created\n")
            dir.create(outdir)
            print(outdir)
        }

        if (is.null(config_file)) {
            cat("No config file provided, creating a default config in", outdir, "\n")
            config_file <- create_config(outdir)
        }

        # argument verificiation
        config <- jsonlite::fromJSON(config_file)

        if (config$isoform_parameters$downsample_ratio > 1 || config$isoform_parameters$downsample_ratio <= 0) {
            stop("downsample_ratio should be between 0 and 1")
        }
        if (!is.null(fastq) &&
            any(!file.exists(fastq))) {
            stop(paste0("Make sure ", fastq, " exists."))
        }
        if (!file.exists(annotation)) {
            stop(paste0("Make sure ", annotation, " exists."))
        }
        if (!file.exists(genome_fa)) {
            stop(paste0("Make sure ", genome_fa, " exists."))
        }

        if (!is.null(genome_bam)) {
            if (any(!file.exists(genome_bam))) {
                stop("Make sure genome_bam exists")
            }
        }

        if (config$pipeline_parameters$do_genome_alignment || config$pipeline_parameters$do_read_realignment) {
            minimap2_dir <- locate_minimap2_dir(minimap2_dir = minimap2_dir)
        }

        if (config$pipeline_parameters$bambu_isoform_identification) {
            if (Matrix::tail(stringr::str_split(annotation, "\\.")[[1]], n = 1) != "gtf") {
                stop("Bambu requires GTF format for annotation file.\n")
            }
        }
        
        n_cores <- parallel::detectCores()
        if (!is.na(n_cores) && config$pipeline_parameters$threads > n_cores) {
                cat("Configured to use", config$pipeline_parameters$threads, "cores, detected", n_cores, "\n")
        }

        return(list(config = config, minimap2_dir = minimap2_dir))
    }
