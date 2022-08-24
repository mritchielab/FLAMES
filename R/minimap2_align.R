#' Minimap2 Align to Genome
#'
#' @description
#' Uses minimap2 to align sequences agains a reference databse.
#' Uses options "-ax splice -t 12 -k14 --secondary=no \code{fa_file} \code{fq_in}"
#'
#' @param config Parsed list of FLAMES config file
#' @param fa_file Path to the fasta file used as a reference database for alignment
#' @param fq_in File path to the fastq file used as a query sequence file
#' @param annot Genome annotation file used to create junction bed files
#' @param outdir Output folder
#' @param minimap2_dir Path to the directory containing minimap2
#' @param prefix String, the prefix (e.g. sample name) for the outputted BAM file
#' @param threads Integer, threads for minimap2 to use, see minimap2 documentation for details,
#' FLAMES will try to detect cores if this parameter is not provided.
#'
#' @return Path to the alignment file
#' @seealso [minimap2_realign()]
#'
#' @importFrom parallel detectCores
#' @importFrom Rsamtools sortBam indexBam asBam
#' @export
#' @examples
#' temp_path <- tempfile()
#' bfc <- BiocFileCache::BiocFileCache(temp_path, ask = FALSE)
#' file_url <- "https://raw.githubusercontent.com/OliverVoogd/FLAMESData/master/data"
#' fastq1 <- bfc[[names(BiocFileCache::bfcadd(bfc, "Fastq1", paste(file_url, "fastq/sample1.fastq.gz", sep = "/")))]]
#' genome_fa <- bfc[[names(BiocFileCache::bfcadd(bfc, "genome.fa", paste(file_url, "SIRV_isoforms_multi-fasta_170612a.fasta", sep = "/")))]]
#' annotation <- bfc[[names(BiocFileCache::bfcadd(bfc, "annot.gtf", paste(file_url, "SIRV_isoforms_multi-fasta-annotation_C_170612a.gtf", sep = "/")))]]
#' outdir <- tempfile()
#' dir.create(outdir)
#' if (is.character(locate_minimap2_dir())) {
#'     minimap2_align(
#'         config = jsonlite::fromJSON(system.file("extdata/SIRV_config_default.json", package = "FLAMES")),
#'         fa_file = genome_fa,
#'         fq_in = fastq1,
#'         annot = annotation,
#'         outdir = outdir
#'     )
#' }
minimap2_align <- function(config, fa_file, fq_in, annot, outdir, minimap2_dir, prefix = NULL, threads = NULL) {
    # if (is.null(threads) && base::system("getconf _NPROCESSORS_ONLN", ignore.stderr = TRUE, ignore.stdout = TRUE) == 0) {
    #    available_cores <- base::system2(command = "getconf", args = c("NPROCESSORS_ONLN"), stdout = TRUE, stderr = TRUE)
    #    if (length(available_cores) == 1 && base::grepl("^[0-9]*$", available_cores)) {
    #        threads <- as.integer(available_cores)
    #    } else {
    #        threads <- 12
    #    }
    # } else {
    #    threads <- 12
    # }
    if (missing("threads") || is.null(threads)) {
        threads <- parallel::detectCores()
    }

    if (!is.null(prefix)) {
        prefix <- paste0(prefix, "_")
    }

    if (missing("minimap2_dir") || is.null(minimap2_dir) || is.na(minimap2_dir)) {
        minimap2_dir <- locate_minimap2_dir()
    }

    samtools <- locate_samtools()

    minimap2_args <- c("-ax", "splice", "-t", threads, "-k14", "--secondary=no", "--seed", config$pipeline_parameters$seed)
    if (config$alignment_parameters$no_flank) {
        minimap2_args <- base::append(minimap2_args, "--splice-flank=no")
    }

    # k8 paftools.js gff2bend gff > bed12
    if (config$alignment_parameters$use_junctions) {
        if (file_test("-f", file.path(minimap2_dir, "k8")) && file_test("-f", file.path(minimap2_dir, "paftools.js"))) {
            paftoolsjs_path <- minimap2_dir
        } else if ((file_test("-f", file.path(dirname(minimap2_dir), "k8")) && file_test("-f", file.path(dirname(minimap2_dir), "paftools.js")))) {
            paftoolsjs_path <- dirname(minimap2_dir)
        } else {
            stop("Could not locate k8 and/or paftools.js in the minimap2 folder, they are required for converting annotation to bed12 files")
        }
        paftoolsjs_status <- base::system2(
            command = file.path(paftoolsjs_path, "k8"),
            args = c(file.path(paftoolsjs_path, "paftools.js"), "gff2bed", annot, ">", file.path(outdir, "tmp_splice_anno.bed12"))
        )
        if (!is.null(base::attr(paftoolsjs_status, "status")) && base::attr(paftoolsjs_status, "status") != 0) {
            stop(paste0("error running k8 paftools.js gff2bed:\n", paftoolsjs_status))
        }
        minimap2_args <- base::append(minimap2_args, c("--junc-bed", file.path(outdir, "tmp_splice_anno.bed12"), "--junc-bonus", "1"))
    }

    # /bin/minimap2 -ax splice -t 12 --junc-bed /.../FLAMES_out/tmp_splice_anno.bed12 --junc-bonus 1 -k14 --secondary=no -o /.../FLAMES_datasets/MuSC/FLAMES_out/tmp_align.sam --seed 2022 /.../GRCm38.primary_assembly.genome.fa /.../trimmed_MSC.fastq.gz
    if (is.character(samtools)) {
        minimap2_status <- base::system2(
            command = file.path(minimap2_dir, "minimap2"),
            args = base::append(minimap2_args, c(
                fa_file, fq_in, "|", samtools, "view -bS -@ 4 -m 2G -o",
                file.path(outdir, paste0(prefix, "tmp_align.bam")), "-"
            ))
        )
        if (!is.null(base::attr(minimap2_status, "status")) && base::attr(minimap2_status, "status") != 0) {
            stop(paste0("error running minimap2:\n", minimap2_status))
        }
        Rsamtools::sortBam(file.path(outdir, paste0(prefix, "tmp_align.bam")), file.path(outdir, paste0(prefix, "align2genome")))
        Rsamtools::indexBam(file.path(outdir, paste0(prefix, "align2genome.bam")))
    } else {
        minimap2_status <- base::system2(
            command = file.path(minimap2_dir, "minimap2"),
            args = base::append(minimap2_args, c(fa_file, fq_in, "-o", file.path(outdir, paste0(prefix, "tmp_align.sam"))))
        )
        if (!is.null(base::attr(minimap2_status, "status")) && base::attr(minimap2_status, "status") != 0) {
            stop(paste0("error running minimap2:\n", minimap2_status))
        }

        Rsamtools::asBam(file.path(outdir, paste0(prefix, "tmp_align.sam")), file.path(outdir, paste0(prefix, "tmp_align")))
        Rsamtools::sortBam(file.path(outdir, paste0(prefix, "tmp_align.bam")), file.path(outdir, paste0(prefix, "align2genome")))
        Rsamtools::indexBam(file.path(outdir, paste0(prefix, "align2genome.bam")))
        file.remove(file.path(outdir, paste0(prefix, "tmp_align.sam")))
        file.remove(file.path(outdir, paste0(prefix, "tmp_align.bam")))
        file.remove(file.path(outdir, paste0(prefix, "tmp_align.bam.bai")))
    }
    if (config$alignment_parameters$use_junctions) {
        file.remove(file.path(outdir, "tmp_splice_anno.bed12"))
    }
    return(file.path(outdir, paste0(prefix, "align2genome.bam")))
}


#' Minimap2 re-align reads to transcriptome
#'
#' @description
#' Uses minimap2 to re-align reads to transcriptome
#'
#' @param config Parsed list of FLAMES config file
#' @param fq_in File path to the fastq file used as a query sequence file
#' @param outdir Output folder
#' @param minimap2_dir Path to the directory containing minimap2
#' @param prefix String, the prefix (e.g. sample name) for the outputted BAM file
#' @param threads Integer, threads for minimap2 to use, see minimap2 documentation for details,
#' FLAMES will try to detect cores if this parameter is not provided.
#'
#' @return Path to the alignment file
#' @seealso [minimap2_align()]
#'
#' @importFrom parallel detectCores
#' @importFrom Rsamtools sortBam indexBam asBam
#' @export
#' @examples
#' temp_path <- tempfile()
#' bfc <- BiocFileCache::BiocFileCache(temp_path, ask = FALSE)
#' file_url <- "https://raw.githubusercontent.com/OliverVoogd/FLAMESData/master/data"
#' fastq1 <- bfc[[names(BiocFileCache::bfcadd(bfc, "Fastq1", paste(file_url, "fastq/sample1.fastq.gz", sep = "/")))]]
#' genome_fa <- bfc[[names(BiocFileCache::bfcadd(bfc, "genome.fa", paste(file_url, "SIRV_isoforms_multi-fasta_170612a.fasta", sep = "/")))]]
#' annotation <- bfc[[names(BiocFileCache::bfcadd(bfc, "annot.gtf", paste(file_url, "SIRV_isoforms_multi-fasta-annotation_C_170612a.gtf", sep = "/")))]]
#' outdir <- tempfile()
#' dir.create(outdir)
#' if (is.character(locate_minimap2_dir())) {
#'     fasta <- annotation_to_fasta(annotation, genome_fa, outdir)
#'     minimap2_realign(
#'         config = jsonlite::fromJSON(system.file("extdata/SIRV_config_default.json", package = "FLAMES")),
#'         fq_in = fastq1,
#'         outdir = outdir
#'     )
#' }
minimap2_realign <- function(config, fq_in, outdir, minimap2_dir, prefix = NULL, threads = NULL) {
    if (missing("threads") || is.null(threads)) {
        threads <- parallel::detectCores()
    }

    if (!is.null(prefix)) {
        prefix <- paste0(prefix, "_")
    }

    samtools <- locate_samtools()

    if (missing("minimap2_dir") || is.null(minimap2_dir) || is.na(minimap2_dir)) {
        minimap2_dir <- locate_minimap2_dir()
    }

    minimap2_args <- c("-ax", "map-ont", "-p", "0.9", "--end-bonus", "10", "-N", "3", "-t", threads, "--seed", config$pipeline_parameters$seed)

    if (is.character(samtools)) {
        minimap2_status <- base::system2(command = file.path(minimap2_dir, "minimap2"), args = base::append(minimap2_args, c(
            file.path(outdir, "transcript_assembly.fa"), fq_in, "|", samtools, "view -bS -@ 4 -m 2G -o",
            file.path(outdir, paste0(prefix, "tmp_align.bam")), "-"
        )))
        if (!is.null(base::attr(minimap2_status, "status")) && base::attr(minimap2_status, "status") != 0) {
            stop(paste0("error running minimap2:\n", minimap2_status))
        }
        Rsamtools::sortBam(file.path(outdir, paste0(prefix, "tmp_align.bam")), file.path(outdir, paste0(prefix, "realign2transcript")))
        Rsamtools::indexBam(file.path(outdir, paste0(prefix, "realign2transcript.bam")))
    } else {
        minimap2_status <- base::system2(
            command = file.path(minimap2_dir, "minimap2"),
            args = base::append(minimap2_args, c(
                file.path(outdir, "transcript_assembly.fa"),
                fq_in,
                "-o",
                file.path(outdir, paste0(prefix, "tmp_align.sam"))
            ))
        )
        if (!is.null(base::attr(minimap2_status, "status")) && base::attr(minimap2_status, "status") != 0) {
            stop(paste0("error running minimap2:\n", minimap2_status))
        }

        Rsamtools::asBam(file.path(outdir, paste0(prefix, "tmp_align.sam")), file.path(outdir, paste0(prefix, "tmp_align")))
        Rsamtools::sortBam(file.path(outdir, paste0(prefix, "tmp_align.bam")), file.path(outdir, paste0(prefix, "realign2transcript")))
        Rsamtools::indexBam(file.path(outdir, paste0(prefix, "realign2transcript.bam")))

        file.remove(file.path(outdir, paste0(prefix, "tmp_align.sam")))
        file.remove(file.path(outdir, paste0(prefix, "tmp_align.bam")))
        file.remove(file.path(outdir, paste0(prefix, "tmp_align.bam.bai")))
    }
    return(file.path(outdir, paste0(prefix, "realign2transcript.bam")))
}

#' locate parent folder of minimap2
#' @description
#' locate minimap2, or validate that minimap2 exists under minimap2_dir
#' @param minimap2_dir (Optional) folder that includes minimap2
#' @return Path to the parent folder of minimap2 if available, or
#' FALSE if minimap2 could not be found.
#' @examples
#' locate_minimap2_dir()
#' @export
locate_minimap2_dir <- function(minimap2_dir = NULL) {
    if (is.null(minimap2_dir)) {
        which_minimap2 <- base::system2(command = "which", args = c("minimap2"), stdout = TRUE, stderr = TRUE)
        if (!is.null(base::attr(which_minimap2, "status")) && base::attr(which_minimap2, "status") != 0) {
            tryCatch(
                {
                    x <- base::system(command = "/usr/bin/modulecmd bash load minimap2", intern = TRUE)
                    minimap2_dir <- dirname(base::system(command = paste0(x, "which minimap2"), intern = TRUE))
                },
                error = function(x) {
                    # stop(paste0("error finding minimap2:\n", which_minimap2))
                    return(FALSE)
                }
            )
        } else {
            minimap2_dir <- dirname(which_minimap2)
        }
    } else if (file_test("-f", minimap2_dir) && grepl("minimap2$", minimap2_dir)) {
        minimap2_dir <- dirname(minimap2_dir)
    } else if (!(file_test("-d", minimap2_dir) && file_test("-f", file.path(minimap2_dir, "minimap2")))) {
        # stop("Error finding minimap2 in minimap2_dir")
        return(FALSE)
    }
    return(minimap2_dir)
}

locate_samtools <- function() {
    which_samtools <- base::system2(command = "which", args = c("samtools"), stdout = TRUE, stderr = TRUE)
    if (!is.null(base::attr(which_samtools, "status")) && base::attr(which_samtools, "status") != 0) {
        tryCatch(
            {
                x <- base::system(command = "/usr/bin/modulecmd bash load samtools", intern = TRUE)
                return(base::system(command = paste0(x, "which samtools"), intern = TRUE))
            },
            error = function(x) {
                return(FALSE)
            }
        )
    } else {
        return("samtools")
    }
}
