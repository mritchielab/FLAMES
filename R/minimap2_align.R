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
#' @return NULL, BAM file is saved to \code{outdir}
#' @seealso [minimap2_realign()]
#'
#' @importFrom parallel detectCores
#' @importFrom Rsamtools sortBam indexBam asBam
#' @export
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
    if (is.null(threads)) {
        threads <- parallel::detectCores()
    }

    if (!is.null(prefix)) {
        prefix <- paste0(prefix, "_")
    }

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
    if (config$alignment_parameters$use_junctions) {
        file.remove(file.path(outdir, "tmp_splice_anno.bed12"))
    }
    return(NULL)
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
#' @return NULL, BAM file is saved to \code{outdir}
#' @seealso [minimap2_align()]
#'
#' @importFrom parallel detectCores
#' @importFrom Rsamtools sortBam indexBam asBam
#' @export
minimap2_realign <- function(config, fq_in, outdir, minimap2_dir, prefix = NULL, threads = NULL) {
    if (is.null(threads)) {
        threads <- parallel::detectCores()
    }

    if (!is.null(prefix)) {
        prefix <- paste0(prefix, "_")
    }

    minimap2_args <- c("-ax", "map-ont", "-p", "0.9", "--end-bonus", "10", "-N", "3", "-t", threads, "--seed", config$pipeline_parameters$seed)
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
    return(NULL)
}

#' Locate minimap2
#'
#' @description Locate the folder containing minimap2, stop if minimap2 is not
#' available in the environment or in the user provided folder. If the user provided
#'  a path to a file, check if it is minimap2 and return its parent folder.
#'
#' @param minimap2_dir User's input for \code{minimap2_dir}
#' @return Path to folder containing minimap2
locate_minimap2_dir <- function(minimap2_dir = NULL) {
    if (is.null(minimap2_dir)) {
        which_minimap2 <- base::system2(command = "which", args = c("minimap2"), stdout = TRUE, stderr = TRUE)
        if (!is.null(base::attr(which_minimap2, "status")) && base::attr(which_minimap2, "status") != 0) {
            stop(paste0("error finding minimap2:\n", which_minimap2))
        } else {
            minimap2_dir <- dirname(which_minimap2)
        }
    } else if (file_test("-f", minimap2_dir) && grepl("minimap2$", minimap2_dir)) {
        minimap2_dir <- dirname(minimap2_dir)
    } else if (!(file_test("-d", minimap2_dir) && file_test("-f", file.path(minimap2_dir, "minimap2")))) {
        stop("Error finding minimap2 in minimap2_dir")
    }
    return(minimap2_dir)
}
