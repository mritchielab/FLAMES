#' Minimap2 Align to Genome
#'
#' @description
#' Uses minimap2 to align sequences agains a reference databse.
#' Uses options '-ax splice -t 12 -k14 --secondary=no \code{fa_file} \code{fq_in}'
#'
#' @param config Parsed list of FLAMES config file
#' @param fa_file Path to the fasta file used as a reference database for alignment
#' @param fq_in File path to the fastq file used as a query sequence file
#' @param annot Genome annotation file used to create junction bed files
#' @param outdir Output folder
#' @param minimap2 Path to minimap2 binary
#' @param k8 Path to the k8 Javascript shell binary
#' @param prefix String, the prefix (e.g. sample name) for the outputted BAM file
#' @param threads Integer, threads for minimap2 to use, see minimap2 documentation for details,
#' FLAMES will try to detect cores if this parameter is not provided.
#' @param samtools path to the samtools binary, required for large datasets since \code{Rsamtools} does not support \code{CSI} indexing
#'
#' @return a \code{data.frame} summarising the reads aligned
#' @seealso [minimap2_realign()]
#'
#' @importFrom Rsamtools sortBam indexBam asBam
#' @export
#' @examples
#' temp_path <- tempfile()
#' bfc <- BiocFileCache::BiocFileCache(temp_path, ask = FALSE)
#' file_url <- 'https://raw.githubusercontent.com/OliverVoogd/FLAMESData/master/data'
#' fastq1 <- bfc[[names(BiocFileCache::bfcadd(bfc, 'Fastq1', paste(file_url, 'fastq/sample1.fastq.gz', sep = '/')))]]
#' genome_fa <- bfc[[names(BiocFileCache::bfcadd(bfc, 'genome.fa', paste(file_url, 'SIRV_isoforms_multi-fasta_170612a.fasta', sep = '/')))]]
#' annotation <- bfc[[names(BiocFileCache::bfcadd(bfc, 'annot.gtf', paste(file_url, 'SIRV_isoforms_multi-fasta-annotation_C_170612a.gtf', sep = '/')))]]
#' outdir <- tempfile()
#' dir.create(outdir)
#' if (!any(is.na(sys_which(c("minimap2", "k8"))))) {
#'     minimap2_align(
#'         config = jsonlite::fromJSON(system.file('extdata/SIRV_config_default.json', package = 'FLAMES')),
#'         fa_file = genome_fa,
#'         fq_in = fastq1,
#'         annot = annotation,
#'         outdir = outdir
#'     )
#' }
minimap2_align <- function(config, fa_file, fq_in, annot, outdir, minimap2 = NA, k8 = NA, samtools = NA,
  prefix = NULL, threads = 1) {
  cat(format(Sys.time(), "%X %a %b %d %Y"), "minimap2_align\n")

  if (!is.null(prefix)) {
    prefix <- paste0(prefix, "_")
  }

  if (missing("minimap2") || is.null(minimap2) || is.na(minimap2) || minimap2 =="") {
    minimap2 <- sys_which("minimap2")
    if (is.na(minimap2)) {
      stop("minimap2 not found, please make sure it is installed and provide its path as the minimap2 argument")
    }
  }

  if (missing("k8") || is.null(k8) || is.na(k8) || k8 =="") {
    k8 <- sys_which("k8")
    if (is.na(k8)) {
      stop("k8 not found, please make sure it is installed and provide its path as the k8 argument")
    }
  }

  if (missing("samtools") || is.null(samtools) || is.na(samtools) || samtools =="") {
    samtools <- sys_which("samtools")
  }

  minimap2_args <- c("-ax", "splice", "-t", threads, "-k14", "--secondary=no",
    "--seed", config$pipeline_parameters$seed)
  if (config$alignment_parameters$no_flank) {
    minimap2_args <- base::append(minimap2_args, "--splice-flank=no")
  }

  # k8 paftools.js gff2bend gff > bed12
  paftoolsjs <- system.file("paftools.js", package = "FLAMES")
  if (config$alignment_parameters$use_junctions) {
    paftoolsjs_status <- base::system2(command = k8,
      args = c(paftoolsjs, "gff2bed", annot,
        ">", file.path(outdir, "tmp_splice_anno.bed12")))
    if (!is.null(base::attr(paftoolsjs_status, "status")) && base::attr(paftoolsjs_status,
      "status") != 0) {
      stop(paste0("error running k8 paftools.js gff2bed:\n", paftoolsjs_status, "\n",
        "Are you using NCBI GFF3? It is not well supported by minimap2's paftools.js"))
      ))
    }
    minimap2_args <- base::append(minimap2_args, c("--junc-bed", file.path(outdir,
      "tmp_splice_anno.bed12"), "--junc-bonus", "1"))
  }

  # /bin/minimap2 -ax splice -t 12 --junc-bed
  # /.../FLAMES_out/tmp_splice_anno.bed12 --junc-bonus 1 -k14 --secondary=no -o
  # /.../FLAMES_datasets/MuSC/FLAMES_out/tmp_align.sam --seed 2022
  # /.../GRCm38.primary_assembly.genome.fa /.../trimmed_MSC.fastq.gz
  if (!is.na(samtools)) {
    minimap2_status <- base::system2(command = minimap2,
      args = base::append(minimap2_args, c(fa_file, fq_in, "|", samtools, "view -bS -@ 4 -o",
        file.path(outdir, paste0(prefix, "tmp_align.bam")), "-")))
    if (!is.null(base::attr(minimap2_status, "status")) && base::attr(minimap2_status,
      "status") != 0) {
      stop(paste0("error running minimap2:\n", minimap2_status))
    }
    sort_status <- base::system2(command = samtools, args = c("sort", file.path(outdir,
      paste0(prefix, "tmp_align.bam")), "-o", file.path(outdir, paste0(prefix,
      "align2genome.bam"))))
    index_status <- base::system2(command = samtools, args = c("index", file.path(outdir,
      paste0(prefix, "align2genome.bam"))))
  } else {
    message("samtools not found, using Rsamtools instead")
    minimap2_status <- base::system2(command = minimap2,
      args = base::append(minimap2_args, c(fa_file, fq_in, "-o", file.path(outdir,
        paste0(prefix, "tmp_align.sam")))))
    if (!is.null(base::attr(minimap2_status, "status")) && base::attr(minimap2_status,
      "status") != 0) {
      stop(paste0("error running minimap2:\n", minimap2_status))
    }

    Rsamtools::asBam(file.path(outdir, paste0(prefix, "tmp_align.sam")), file.path(outdir,
      paste0(prefix, "tmp_align")))
    Rsamtools::sortBam(file.path(outdir, paste0(prefix, "tmp_align.bam")), file.path(outdir,
      paste0(prefix, "align2genome")))
    Rsamtools::indexBam(file.path(outdir, paste0(prefix, "align2genome.bam")))
    file.remove(file.path(outdir, paste0(prefix, "tmp_align.sam")))
  }
  file.remove(file.path(outdir, paste0(prefix, "tmp_align.bam")))

  if (config$alignment_parameters$use_junctions) {
    file.remove(file.path(outdir, "tmp_splice_anno.bed12"))
  }
  return(get_flagstat(file.path(outdir, paste0(prefix, "align2genome.bam")), samtools))
}


#' Minimap2 re-align reads to transcriptome
#'
#' @description
#' Uses minimap2 to re-align reads to transcriptome
#'
#' @param config Parsed list of FLAMES config file
#' @param fq_in File path to the fastq file used as a query sequence file
#' @param outdir Output folder
#' @param minimap2 Path to minimap2 binary
#' @param prefix String, the prefix (e.g. sample name) for the outputted BAM file
#' @param threads Integer, threads for minimap2 to use, see minimap2 documentation for details,
#' FLAMES will try to detect cores if this parameter is not provided.
#'
#' @return a \code{data.frame} summarising the reads aligned
#' @seealso [minimap2_align()]
#'
#' @importFrom Rsamtools sortBam indexBam asBam
#' @export
#' @examples
#' temp_path <- tempfile()
#' bfc <- BiocFileCache::BiocFileCache(temp_path, ask = FALSE)
#' file_url <- 'https://raw.githubusercontent.com/OliverVoogd/FLAMESData/master/data'
#' fastq1 <- bfc[[names(BiocFileCache::bfcadd(bfc, 'Fastq1', paste(file_url, 'fastq/sample1.fastq.gz', sep = '/')))]]
#' genome_fa <- bfc[[names(BiocFileCache::bfcadd(bfc, 'genome.fa', paste(file_url, 'SIRV_isoforms_multi-fasta_170612a.fasta', sep = '/')))]]
#' annotation <- bfc[[names(BiocFileCache::bfcadd(bfc, 'annot.gtf', paste(file_url, 'SIRV_isoforms_multi-fasta-annotation_C_170612a.gtf', sep = '/')))]]
#' outdir <- tempfile()
#' dir.create(outdir)
#' if (!any(is.na(sys_which(c("minimap2", "k8"))))) {
#'     fasta <- annotation_to_fasta(annotation, genome_fa, outdir)
#'     minimap2_realign(
#'         config = jsonlite::fromJSON(system.file('extdata/SIRV_config_default.json', package = 'FLAMES')),
#'         fq_in = fastq1,
#'         outdir = outdir
#'     )
#' }
minimap2_realign <- function(config, fq_in, outdir, minimap2, samtools = NULL, prefix = NULL,
  threads = 1) {
  cat(format(Sys.time(), "%X %a %b %d %Y"), "minimap2_realign\n")

  if (!is.null(prefix)) {
    prefix <- paste0(prefix, "_")
  }

  if (missing("minimap2") || is.null(minimap2) || is.na(minimap2) || minimap2 =="") {
    minimap2 <- sys_which("minimap2")
    if (is.na(minimap2)) {
      stop("minimap2 not found, please make sure it is installed and provide its path as the minimap2 argument")
    }
  }

  if (missing("samtools") || is.null(samtools) || is.na(samtools) || samtools =="") {
    samtools <- sys_which("samtools")
  }

  minimap2_args <- c("-ax", "map-ont", "-p", "0.9", "--end-bonus", "10", "-N",
    "3", "-t", threads, "--seed", config$pipeline_parameters$seed)

  if (!is.na(samtools)) {
    minimap2_status <- base::system2(command = minimap2,
      args = base::append(minimap2_args, c(file.path(outdir, "transcript_assembly.fa"),
        fq_in, "|", samtools, "view -bS -@ 4 -o", file.path(outdir,
          paste0(prefix, "tmp_align.bam")), "-")))
    if (!is.null(base::attr(minimap2_status, "status")) && base::attr(minimap2_status,
      "status") != 0) {
      stop(paste0("error running minimap2:\n", minimap2_status))
    }
    sort_status <- base::system2(command = samtools, args = c("sort", file.path(outdir,
      paste0(prefix, "tmp_align.bam")), "-o", file.path(outdir, paste0(prefix,
      "realign2transcript.bam"))))
    index_status <- base::system2(command = samtools, args = c("index", file.path(outdir,
      paste0(prefix, "realign2transcript.bam"))))
  } else {
    message("samtools not found, using Rsamtools instead")
    minimap2_status <- base::system2(command = minimap2,
      args = base::append(minimap2_args, c(file.path(outdir, "transcript_assembly.fa"),
        fq_in, "-o", file.path(outdir, paste0(prefix, "tmp_align.sam")))))
    if (!is.null(base::attr(minimap2_status, "status")) && base::attr(minimap2_status,
      "status") != 0) {
      stop(paste0("error running minimap2:\n", minimap2_status))
    }

    Rsamtools::asBam(file.path(outdir, paste0(prefix, "tmp_align.sam")), file.path(outdir,
      paste0(prefix, "tmp_align")))
    Rsamtools::sortBam(file.path(outdir, paste0(prefix, "tmp_align.bam")), file.path(outdir,
      paste0(prefix, "realign2transcript")))
    Rsamtools::indexBam(file.path(outdir, paste0(prefix, "realign2transcript.bam")))

    file.remove(file.path(outdir, paste0(prefix, "tmp_align.sam")))
  }
  file.remove(file.path(outdir, paste0(prefix, "tmp_align.bam")))

  return(get_flagstat(file.path(outdir, paste0(prefix, "realign2transcript.bam")),
    samtools))
}

#' Sys.which wrapper
#' Wrapper for Sys.which that replaces "" with NA
#' @description
#' The \code{base::Sys.which} function returns "" if the command is not found
#' on some systems and \code{NA} on others. This wrapper replaces "" with
#' \code{NA}
#' @param command character, the command to search for
#' @return character, the path to the command or \code{NA}
#' @examples
#' sys_which("minimap2")
#' @export
sys_which <- function(command) {
  which_command <- Sys.which(command)
  # replace "" with NA
  which_command[which_command == ""] <- NA
  return(which_command)
}

# total mapped primary secondary
get_flagstat <- function(bam, samtools_path) {
  stats_df <- data.frame(total = 0, mapped = 0, primary = 0, secondary = 0)
  rownames(stats_df) <- bam
  if (!missing("samtools_path") && is.character(samtools_path)) {
    output <- base::system2(samtools_path, c("flagstat", bam), stdout = TRUE)
    stats_df["total"] <- as.numeric(regmatches(output[grepl("in total", output)],
      regexpr("\\d+", output[grepl("in total", output)]))[1])
    stats_df["mapped"] <- as.numeric(regmatches(output[grepl("mapped", output)],
      regexpr("(?<!S)\\d+", output[grepl("mapped", output)], perl = TRUE))[1])
    stats_df["primary"] <- as.numeric(regmatches(output[grepl("primary$", output)],
      regexpr("(?<!S)\\d+", output[grepl("primary$", output)], perl = TRUE))[1])
    stats_df["secondary"] <- as.numeric(regmatches(output[grepl("secondary$",
      output)], regexpr("(?<!S)\\d+", output[grepl("secondary$", output)],
      perl = TRUE))[1])
  } else {
    output <- utils::capture.output(Rsamtools::quickBamFlagSummary(bam))
    stats_df["total"] <- as.numeric(regmatches(output[grepl("^All records", output)],
      regexpr("\\d+", output[grepl("^All records", output)]))[1])
    stats_df["mapped"] <- as.numeric(regmatches(output[grepl("record is mapped",
      output)], regexpr("(?<!S)\\d+", output[grepl("record is mapped", output)],
      perl = TRUE))[1])
    stats_df["primary"] <- as.numeric(regmatches(output[grepl("primary alignment",
      output)], regexpr("(?<!S)\\d+", output[grepl("primary alignment", output)],
      perl = TRUE))[1])
    stats_df["secondary"] <- as.numeric(regmatches(output[grepl("secondary alignment",
      output)], regexpr("(?<!S)\\d+", output[grepl("primary alignment", output)],
      perl = TRUE))[1])
  }
  stats_df
}

# generic pie chart
#' @importFrom ggplot2 ggplot aes geom_bar ggtitle coord_polar element_blank theme position_stack theme_bw geom_histogram ggtitle ylab xlab geom_text
plot_flagstat <- function(flagstat) {
  flagstat[["unmapped"]] <- flagstat[["total"]] - flagstat[["mapped"]]
  tidyr::pivot_longer(as.data.frame(flagstat), everything()) %>%
    dplyr::filter(!name %in% c("total", "mapped")) %>%
    ggplot(aes(x = "", y = value, label = value, fill = name)) + geom_bar(stat = "identity") +
    coord_polar("y") + ggtitle("Alignment summary") + geom_text(position = position_stack(vjust = 0.5)) +
    labs(x = NULL, y = NULL) + theme_bw() + theme(panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_blank(),
    axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(),
    axis.ticks.y = element_blank())
}
