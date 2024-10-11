#' filter annotation for plotting coverages
#'
#' @description Removes isoform annotations that could produce ambigious reads, such as isoforms
#' that only differ by the 5' / 3' end. This could be useful for plotting average
#' coverage plots.
#'
#' @importFrom txdbmaker makeTxDbFromGFF makeTxDbFromGRanges
#' @importFrom rtracklayer import
#' @importFrom S4Vectors split
#' @importFrom GenomicRanges strand
#' @importFrom BiocGenerics start end
#'
#' @param annotation path to the GTF annotation file, or the parsed GenomicRanges
#' object.
#' @param keep string, one of 'tss_differ' (only keep isoforms that all differ
#' by the transcription start site position), 'tes_differ' (only keep those that
#' differ by the transcription end site position), 'both' (only keep those that
#' differ by both the start and end site), or 'single_transcripts' (only keep
#' genes that contains a sinlge transcript).
#' @return GenomicRanges of the filtered isoforms
#' @examples
#' filtered_annotation <- filter_annotation(
#'   system.file("extdata", "rps24.gtf.gz", package = 'FLAMES'), keep = 'tes_differ')
#' filtered_annotation
#'
#' @md
#' @export
filter_annotation <- function(annotation, keep = "tss_differ") {
  if (is.character(annotation)) {
    annotation <- annotation |>
      txdbmaker::makeTxDbFromGFF() |>
      GenomicFeatures::transcripts()
  } else {
    annotation <- annotation |>
      txdbmaker::makeTxDbFromGRanges() |>
      GenomicFeatures::transcripts()
  }

  unique_fn <- function(x, keep) {
    if (keep == "tss_differ") {
      return(!duplicated(GenomicRanges::start(x)) & !duplicated(GenomicRanges::start(x), fromLast = TRUE))
    }
    if (keep == "tes_differ") {
      return(!duplicated(GenomicRanges::end(x)) & !duplicated(GenomicRanges::end(x), fromLast = TRUE))
    }
    if (keep == "both") {
      return(unique_fn(x, "tss_differ") & unique_fn(x, "tes_differ"))
    }
  }

  return(annotation[unique_fn(annotation, keep)])
}

#' plot read coverages
#'
#' @description Plot the average read coverages for each length bin or a
#' perticular isoform
#'
#' @importFrom txdbmaker makeTxDbFromGFF
#' @importFrom Rsamtools ScanBamParam
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr filter mutate group_by summarize_at summarise across arrange left_join
#' @importFrom ggplot2 ggplot geom_line aes xlab ylab guides theme_minimal guide_legend
#' @importFrom stats weighted.mean
#'
#' @param bam, path to the BAM file (aligning reads to the transcriptome), or
#' the (GenomicAlignments::readGAlignments) parsed GAlignments object
#' @param isoform string vector, provide isoform names to plot the coverage for the
#' corresponding isoforms, or provide NULL to plot average coverages for each
#' length bin
#' @param quantiles numeric vector to specify the quantiles to bin the transcripts lengths
#' by if length_bins is missing. The length bins will be determined such that the read
#' counts are distributed acording to the quantiles.
#' @param length_bins, numeric vector to specify the sizes to bin the isoforms by
#' @param weight_fn "read_counts" or "sigmoid", determins how the transcripts
#' should be weighted within length bins.
#' @param annotation path to the GTF annotation file, or the parsed txdb object
#' @param coding_only logical, if TRUE, only plot the coverage for coding transcripts
#' @param remove_utr logical, if TRUE, trim the UTRs from the transcripts
#' @return a ggplot2 object of the coverage plot(s)
#' @examples
#' temp_path <- tempfile()
#' bfc <- BiocFileCache::BiocFileCache(temp_path, ask = FALSE)
#' file_url <- 'https://raw.githubusercontent.com/OliverVoogd/FLAMESData/master/data'
#' fastq1 <- bfc[[names(BiocFileCache::bfcadd(bfc, 'Fastq1', paste(file_url, 'fastq/sample1.fastq.gz', sep = '/')))]]
#' genome_fa <- bfc[[names(BiocFileCache::bfcadd(bfc, 'genome.fa', paste(file_url, 'SIRV_isoforms_multi-fasta_170612a.fasta', sep = '/')))]]
#' annotation <- bfc[[names(BiocFileCache::bfcadd(bfc, 'annot.gtf', paste(file_url, 'SIRV_isoforms_multi-fasta-annotation_C_170612a.gtf', sep = '/')))]]
#' outdir <- tempfile()
#' dir.create(outdir)
#' fasta <- annotation_to_fasta(annotation, genome_fa, outdir)
#' minimap2_realign(
#'   config = jsonlite::fromJSON(
#'     system.file("extdata", "config_sclr_nanopore_3end.json", package = "FLAMES")),
#'   fq_in = fastq1,
#'   outdir = outdir
#' )
#' plot_coverage(bam = file.path(outdir, 'realign2transcript.bam'))
#' @md
#' @export
plot_coverage <- function(bam, isoform = NULL, quantiles = c(0, 0.2375, 0.475, 0.7125, 0.95, 1),
    length_bins = c(0, 1, 2, 5, 10, Inf), weight_fn = "read_counts", annotation, coding_only = FALSE, remove_utr = FALSE) {

  if (coding_only || remove_utr) {
    if (missing(annotation)) {
      stop("Please provide the annotation file to remove UTRs or select for coding transcripts.")
    } else if (is.character(annotation)) {
      annotation <- txdbmaker::makeTxDbFromGFF(annotation)
      transcript_info <- transcript_coverage(bam, isoform, weight_fn, annotation, coding_only, remove_utr)
    } else {
      transcript_info <- transcript_coverage(bam, isoform, weight_fn, annotation, coding_only, remove_utr)
    }
  } else {
    transcript_info <- transcript_coverage(bam, isoform, weight_fn)
  }

  if (!is.null(isoform)) {
    p <- transcript_info |>
      tidyr::pivot_longer(paste0("coverage_", 1:100), names_to = "x", values_to = "coverage") |>
      dplyr::mutate(x = as.numeric(gsub("coverage_", "", x))) |>
      ggplot2::ggplot(aes(x = x, y = coverage, color = transcript)) + 
      geom_line() +
      ggplot2::xlab("% Position in transcript") +
      ggplot2::ylab("Coverage") +
      ggplot2::guides(color = ggplot2::guide_legend(title = "Transcript"))
    return(p)
  }

  if (missing(length_bins)) {
    message("Using quantiles to bin transcripts.")
    # length_bins <- c(0, quantile(transcript_info$tr_length, quantiles), Inf) / 1000
    # message("Length bins: ", paste0(length_bins, collapse = ", "))
    bin_ranges <- transcript_info |>
      dplyr::group_by(tr_length) |>
      dplyr::summarise(total_counts = sum(read_counts), .groups = "drop") |>
      dplyr::arrange(tr_length) |>
      dplyr::mutate(cumpct = cumsum(total_counts) / sum(total_counts)) |>
      dplyr::mutate(length_bin = cut(cumpct, breaks = quantiles, include.lowest = TRUE))

    transcript_info <- transcript_info |>
      dplyr::left_join(bin_ranges, by = "tr_length", relationship = "many-to-one") |>
      dplyr::group_by(length_bin) |>
      dplyr::mutate(
        min_length = round(min(tr_length)/1000, 2),
        max_length = round(max(tr_length)/1000, 2),
        length_bin = paste0(min_length, '-', max_length, ' kb, ', length_bin)
      ) |>
      dplyr::ungroup()
  } else {
    message("Using fixed length bins.")
    transcript_info$length_bin <- cut(transcript_info$tr_length / 1000, breaks = length_bins)
    levels(transcript_info$length_bin) <- 
      paste0(length_bins[-length(length_bins)], '-', length_bins[-1], ' kb')
  }

  mean_coverage <- transcript_info |>
    dplyr::group_by(length_bin) |>
    dplyr::summarise(dplyr::across(paste0("coverage_", 1:100), ~ stats::weighted.mean(.,
      w = weight)))

  p <- mean_coverage |>
    tidyr::pivot_longer(paste0("coverage_", 1:100), names_to = "x", values_to = "coverage") |>
    dplyr::mutate(x = as.numeric(gsub("coverage_", "", x))) |>
    ggplot2::ggplot(aes(x = x, y = coverage, color = length_bin)) +
    ggplot2::geom_line() +
    ggplot2::xlab("% Position in transcript") +
    ggplot2::ylab("Weighted average coverage") +
    ggplot2::theme_minimal()
  if (missing(length_bins)) {
    p <- p + ggplot2::guides(color = ggplot2::guide_legend(title = "Transcript length, read count quantile"))
  } else {
    p <- p + ggplot2::guides(color = ggplot2::guide_legend(title = "Transcript length bins"))
  }

  return(p)
}


#' @importFrom tibble as_tibble
#' @importFrom GenomeInfoDb seqlengths
#' @importFrom Rsamtools ScanBamParam
#' @importFrom GenomicFeatures cdsBy threeUTRsByTranscript fiveUTRsByTranscript
#' @importFrom GenomicAlignments readGAlignments seqnames
#' @importFrom GenomicRanges width strand granges coverage
transcript_coverage <- function(bam, isoform = NULL, weight_fn = "read_counts", annotation, coding_only = FALSE, remove_utr = FALSE) {

  if (!is(bam, "GAlignments")) {
    bam <- GenomicAlignments::readGAlignments(bam, param = Rsamtools::ScanBamParam(mapqFilter = 5))
  }

  if (!is.null(isoform)) {
    bam <- bam[GenomicAlignments::seqnames(bam) %in% isoform]
  }

  read_counts <- table(GenomicAlignments::seqnames(bam))
  transcript_names <- names(read_counts)

  cover <- bam |>
    GenomicRanges::granges() |>
    GenomicRanges::coverage()
  cover <- cover[transcript_names]
  
  if (coding_only) {
    message("Filtering for coding transcripts...")
    transcript_names <- intersect(
      transcript_names,
      names(GenomicFeatures::cdsBy(annotation, by="tx", use.names=TRUE))
    )
    cover <- cover[transcript_names]
  }

  if (remove_utr) {
    message("Getting UTRs...")
    threeUTR.data <- GenomicFeatures::threeUTRsByTranscript(annotation, use.names=T)
    fiveUTR.data <- GenomicFeatures::fiveUTRsByTranscript(annotation, use.names=T)
    threeUTR.data <- threeUTR.data[names(threeUTR.data) %in% transcript_names]
    fiveUTR.data <- fiveUTR.data[names(fiveUTR.data) %in% transcript_names]

    message("Calculating UTR length...")
    threeUTR.len <- rep(0, length(transcript_names))
    fiveUTR.len <- rep(0, length(transcript_names))
    names(threeUTR.len) <- transcript_names
    names(fiveUTR.len) <- transcript_names
    threeUTR.len[names(threeUTR.data)] <- threeUTR.data |> width() |> sum()
    fiveUTR.len[names(fiveUTR.data)] <- fiveUTR.data |> width() |> sum()

    message("Removing UTRs...")
    cover <-
      lapply(seq_along(cover), function(x) {
        cover[[x]][(fiveUTR.len[x] + 1): (length(cover[[x]])-threeUTR.len[x])]
      })
    names(cover) <- transcript_names

    transcript_info <- data.frame(
      tr_length = GenomeInfoDb::seqlengths(bam)[transcript_names] - threeUTR.len - fiveUTR.len,
      read_counts = as.data.frame(read_counts[transcript_names])$Freq
    )
  } else {
    transcript_info <- data.frame(tr_length = GenomeInfoDb::seqlengths(bam)[transcript_names],
      read_counts = as.data.frame(read_counts[transcript_names])$Freq)
  }

  cover <- cover |>
    sapply(function(x) {
      x[round(seq(1, length(x), length.out = 100), 0)] |>
        as.integer()
    }) |>
    t() |>
    as.data.frame()
  colnames(cover) <- paste0("coverage_", 1:100)

  if (weight_fn == "sigmoid") {
    weight_fn <- function(mat, read_counts) {
      sigmoid <- function(x) {
        exp(x) / (exp(x) + 1)
      }
      sigmoid((read_counts - 2000) / 500)
    }
  } else if (weight_fn == "read_counts") {
    weight_fn <- function(mat, read_counts) {read_counts}
  }

  cover <- cover[transcript_names,]
  cover <- cover / transcript_info$read_counts # scale by read counts
  transcript_info <- cbind(transcript_info, cover)
  transcript_info$weight <- weight_fn(mat, transcript_info$read_counts)
  transcript_info <- transcript_info |>
    tibble::as_tibble(rownames = "transcript")
  return(transcript_info)
}
