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
#' @importFrom future plan
wrt_tr_to_csv <-
    function(bc_tr_count_dict,
             transcript_dict,
             csv_f,
             transcript_dict_ref = NULL,
             has_UMI = TRUE) {
        future::plan(future::multisession)
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

#' Transcript quantification
#' @description Calculate the transcript count matrix by parsing the re-alignment file.
#' @param annotation The file path to the annotation file in GFF3 format
#' @param outdir The path to directory to store all output files.
#' @param config Parsed FLAMES configurations.
#' @param pipeline The pipeline type as a character string, either \code{sc_single_sample} (single-cell, single-sample),
#' \code{bulk} (bulk, single or multi-sample), or \code{sc_multi_sample} (single-cell, multiple samples)
#' @return The count matrix will be saved in the output folder as \code{transcript_count.csv.gz}.
#' @importFrom reticulate import_from_path dict
#' @examples
#' temp_path <- tempfile()
#' bfc <- BiocFileCache::BiocFileCache(temp_path, ask = FALSE)
#' file_url <- "https://raw.githubusercontent.com/OliverVoogd/FLAMESData/master/data"
#' fastq1 <- bfc[[names(BiocFileCache::bfcadd(bfc, "Fastq1", paste(file_url, "fastq/sample1.fastq.gz", sep = "/")))]]
#' genome_fa <- bfc[[names(BiocFileCache::bfcadd(bfc, "genome.fa", paste(file_url, "SIRV_isoforms_multi-fasta_170612a.fasta", sep = "/")))]]
#' annotation <- bfc[[names(BiocFileCache::bfcadd(bfc, "annot.gtf", paste(file_url, "SIRV_isoforms_multi-fasta-annotation_C_170612a.gtf", sep = "/")))]]
#' outdir <- tempfile()
#' dir.create(outdir)
#' fasta <- annotation_to_fasta(annotation, genome_fa, outdir)
#' config <- jsonlite::fromJSON(create_config(outdir, bambu_isoform_identification = TRUE, min_tr_coverage = 0.1, min_read_coverage = 0.1, min_sup_cnt = 1))
#' file.copy(annotation, file.path(outdir, "isoform_annotated.gtf"))
#' \dontrun{
#' if (is.character(locate_minimap2_dir())) {
#'     minimap2_realign(
#'         config = config, outdir = outdir,
#'         fq_in = fastq1
#'     )
#'     quantify(annotation, outdir, config, pipeline = "bulk")
#' }
#' }
#' @export
quantify <- function(annotation, outdir, config, pipeline = "sc_single_sample") {
    cat(format(Sys.time(), "%X %a %b %d %Y"), "quantify\n")

    if (grepl("\\.gff3?(\\.gz)?$", annotation)) {
        warning("Annotation in GFF format may cause errors. Please consider using GTF formats.\n")
    }

    realign_bam <- list.files(outdir)[grepl("_?realign2transcript\\.bam$", list.files(outdir))]
    cat("Found realignment file(s): ")
    cat(paste0("\t", paste(realign_bam, collapse = "\n\t"), "\n"))

    if (length(realign_bam) != 1 && grepl("single_sample", pipeline)) {
        stop("Incorrect number of realignment files found.\n")
    }
    
    callBasilisk(flames_env, function(config_dict, annotation, outdir, pipeline) {
        python_path <- system.file("python", package = "FLAMES")
        count <- reticulate::import_from_path("count_tr", python_path)
        count$quantification(config_dict, annotation, outdir, pipeline)
    },
    config_dict = reticulate::dict(config),
    annotation = annotation,
    outdir = outdir,
    pipeline = pipeline
    )


}

# example for Rsamtools
#' @importFrom Rsamtools BamFile  scanBam isIncomplete
#quantify_tmp <- function(bamFileName) {
#    bf <- Rsamtools::BamFile(bamFileName, yieldSize=100)
#    while (Rsamtools::isIncomplete(bf)) {
#        print(Rsamtools::scanBam(bf)[[1]][[1]])
#        Sys.sleep(3)
#    }
#}
