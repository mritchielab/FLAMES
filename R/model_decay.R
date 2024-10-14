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

read_counts_sigmoid_weight <- function(counts, inflection_idx = 10, inflection_max = 1000, steepness = 5) {
  # sigmoid function
  sigmoid_weight <- function(x, infl, stp) {
    1 / (1 + exp(-stp * (x - infl)))
  }

  if (length(counts) <= inflection_idx) {
    message("The number of transcripts is less than the inflection index, returning equal weights for the current bin.")
    return(rep(1, length(counts)))
  }
  inflection <- sort(counts, decreasing = TRUE)[inflection_idx]
  if (inflection > inflection_max) {
    inflection <- inflection_max
  }
  steepness <- steepness / inflection
  return(sigmoid_weight(counts, inflection, steepness))
}

#' Weight transcripts by read counts
#'
#' @description Given a vector of read counts, return a vector of weights.
#' The weights could be either the read counts themselves (\code{type = 'counts'}),
#' a binary vector of 0s and 1s where 1s are assigned to transcripts with read
#' counts above a threshold (\code{type = 'equal', min_counts = 1000}), or a
#' sigmoid function of the read counts (\code{type = 'sigmoid'}). The sigmoid
#' function is defined as \code{1 / (1 + exp(-steepness/inflection * (x - inflection)))}.
#'
#' @param counts numeric vector of read counts
#' @param type string, one of 'counts', 'sigmoid', or 'equal'
#' @param min_counts numeric, the threshold for the 'equal' type
#' @param inflection_idx numeric, the index of the read counts to determine the
#' inflection point for the sigmoid function. The default is 10, i.e. the 10th
#' highest read count will be the inflection point.
#' @param inflection_max numeric, the maximum value for the inflection point.
#' If the inflection point according to the inflection_idx is higher than this
#' value, the inflection point will be set to this value instead.
#' @param steepness numeric, the steepness of the sigmoid function

#' @return numeric vector of weights
#' @examples
#' weight_transcripts(1:2000)
#' par(mfrow = c(2, 2))
#' plot(
#'   1:2000, weight_transcripts(1:2000, type = 'sigmoid'),
#'   type = 'l', xlab = 'Read counts', ylab = 'Sigmoid weight'
#' )
#' plot(
#'   1:2000, weight_transcripts(1:2000, type = 'counts'),
#'   type = 'l', xlab = 'Read counts', ylab = 'Weight by counts'
#' )
#' plot(
#'   1:2000, weight_transcripts(1:2000, type = 'equal'),
#'   type = 'l', xlab = 'Read counts', ylab = 'Equal weights'
#' )
#'
#' @md
#' @export
weight_transcripts <- function(counts, type = 'sigmoid', min_counts = 1000,
    inflection_idx = 10, inflection_max = 1000, steepness = 5) {
  if (type == 'counts') {
    return(counts)
  } else if (type == 'equal') {
    weights <- ifelse(counts > min_counts, 1, 0)
    return(weights)
  } else if (type == 'sigmoid') {
    return(read_counts_sigmoid_weight(counts, inflection_idx, inflection_max, steepness))
  }
}

#' Convolution filter for smoothing transcript coverages
#'
#' @description Filter out transcripts with sharp drops / rises in coverage,
#' to be used in \code{filter_coverage} to remove transcripts with potential
#' misalignments / internal priming etc. Filtering is done by convolving the
#' coverage with a kernal of 1s and -1s (e.g. \code{c(1, 1, -1, -1), where
#' the width of the 1s and -1s are determined by the \code{width} parameter),
#' and check if the maximum absolute value of the convolution is below a
#' threshold. If the convolution is below the threshold, \code{TRUE} is returned,
#' otherwise \code{FALSE}.
#'
#' @param x numeric vector of coverage values
#' @param threshold numeric, the threshold for the maximum absolute value of the
#' convolution
#' @param width numeric, the width of the 1s and -1s in the kernal. E.g. \code{width = 2}
#' will result in a kernal of \code{c(1, 1, -1, -1)}
#' @param trim numeric, the proportion of the coverage values to ignore at
#' both ends before convolution.
#' @return logical, \code{TRUE} if the transcript passes the filter, \code{FALSE} otherwise
#' @examples
#' # A >30% drop in coverage will fail the filter with threshold = 0.3
#' convolution_filter(c(1, 1, 1, 0.69, 0.69, 0.69), threshold = 0.3)
#' convolution_filter(c(1, 1, 1, 0.71, 0.7, 0.7), threshold = 0.3)
#'
#' @md
#' @export
convolution_filter <- function(x, threshold = 0.15, width = 2, trim = 0.05) {
  # trim off both ends
  threshold <- threshold * width
  trimmed <- x[round(length(x) * trim):round(length(x) * (1 - trim))]
  kernal <- c(rep(1, width), rep(-1, width))
  conv <- convolve(trimmed, kernal, type = "filter")
  return(max(abs(conv)) < threshold)
}

#' Get read coverages from BAM file
#'
#' @description Get the read coverages for each transcript in the BAM file (or a GAlignments object).
#' The read coverages are sampled at 100 positions along the transcript, and the coverage is scaled
#' by dividing the coverage at each position by the total read counts for the transcript. If a BAM
#' file is provided, alignment with MAPQ < 5, secondary alignments and supplementary alignments
#' are filtered out. A \code{GAlignments} object can also be provided in case alternative filtering
#' is desired.
#'
#' @importFrom GenomicAlignments readGAlignments seqnames
#' @importFrom Rsamtools scanBamFlag ScanBamParam
#' @importFrom dplyr filter
#' @param bam path to the BAM file, or a parsed GAlignments object
#' @param remove_UTR logical, if \code{TRUE}, remove the UTRs from the coverage
#' @param annotation (Required if \code{remove_UTR = TRUE}) path to the GTF annotation file
#' @return a tibble of the transcript information and coverages, with the following columns:
#' \itemize{
#' \item \code{transcript}: the transcript name / ID
#' \item \code{read_counts}: the total number of aligments for the transcript
#' \item \code{coverage_1-100}: the coverage at each of the 100 positions along the transcript
#' \item \code{tr_length}: the length of the transcript
#' }
#' @examples
#' # Create a BAM file with minimap2_realign
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
#' x <- get_coverage(file.path(outdir, 'realign2transcript.bam'))
#' head(x)
#' @md
#' @export
get_coverage <- function(bam, remove_UTR = FALSE, annotation) {

  if (!is(bam, "GAlignments")) {
    bam <- GenomicAlignments::readGAlignments(bam,
      param = Rsamtools::ScanBamParam(
        Rsamtools::scanBamFlag(isSecondaryAlignment = FALSE, isSupplementaryAlignment = FALSE),
        mapqFilter = 5
      )
    )
  }

  read_counts_df <- bam |>
    GenomicAlignments::seqnames() |>
    table() |>
    as.data.frame() |>
    dplyr::filter(Freq > 0)
  coverage_rlel <- get_coverage_RleList(bam)
  if (remove_UTR) {
    coverage_rlel <- remove_UTR_coverage(coverage_rlel, annotation)
  }
  transcript_info <- sample_coverage(coverage_rlel, read_counts_df)

  return(transcript_info)
}

#' Filter transcript coverage
#'
#' @description Filter the transcript coverage by applying a filter function to the
#' coverage values.
#'
#' @importFrom dplyr rowwise mutate c_across ungroup pull
#' @param x The tibble returned by \code{\link{get_coverage}}, or a BAM file path, or
#' a GAlignments object.
#' @param filter_fn The filter function to apply to the coverage values. The function
#' should take a numeric vector of coverage values and return a logical value (TRUE if
#' the transcript passes the filter, FALSE otherwise). The default filter function is
#' \code{\link{convolution_filter}}, which filters out transcripts with sharp drops /
#' rises in coverage.
#' @return a tibble of the transcript information and coverages, with transcipts that
#' pass the filter
#' @examples
#' # Create a BAM file with minimap2_realign
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
#' x <- get_coverage(file.path(outdir, 'realign2transcript.bam'))
#' nrow(x)
#' filter_coverage(x) |>
#'   nrow()
#' @md
#' @export
filter_coverage <- function(x, filter_fn = convolution_filter) {
  if (is.character(x) || is(x, "GAlignments")) {
    transcript_info <- get_coverage(x)
  } else {
    transcript_info <- x
  }

  filter_idx <- transcript_info |>
    dplyr::rowwise() |>
    dplyr::mutate(filter_res = filter_fn(dplyr::c_across(starts_with("coverage_")))) |>
    dplyr::ungroup() |>
    dplyr::pull(filter_res)
  total_reads <- transcript_info$read_counts |> sum()
  removed_reads <- transcript_info$read_counts[!filter_idx] |> sum()
  message(length(filter_idx), " transcripts found in the BAM file.")
  message(sum(!filter_idx), "(", round(sum(!filter_idx) / length(filter_idx) * 100, 2),
    "%) transcripts failed the filter.")
  message("Failed transcripts account for ", removed_reads, " reads, out of ",
    total_reads, "(", round(removed_reads / total_reads * 100, 2), "%) reads in total.")

  transcript_info <- transcript_info[filter_idx, ]
  return(transcript_info)
}

#' plot read coverages
#'
#' @description Plot the average read coverages for each length bin or a
#' perticular isoform
#'
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr group_by summarise arrange mutate left_join ungroup across filter
#' @importFrom ggplot2 ggplot geom_line aes xlab ylab guides theme_minimal guide_legend
#' @importFrom stats weighted.mean
#'
#' @param x, path to the BAM file (aligning reads to the transcriptome), or
#' the (GenomicAlignments::readGAlignments) parsed GAlignments object, or the
#' tibble returned by \code{\link{get_coverage}}, or the filtered tibble returned
#' by \code{\link{filter_coverage}}.
#' @param quantiles numeric vector to specify the quantiles to bin the transcripts lengths
#' by if length_bins is missing. The length bins will be determined such that the read
#' counts are distributed acording to the quantiles.
#' @param length_bins, numeric vector to specify the sizes to bin the transcripts by
#' @param weight_fn function to calculate the weights for the transcripts. The function
#' should take a numeric vector of read counts and return a numeric vector of weights.
#' The default function is \code{\link{weight_transcripts}}, you can change its default
#' parameters by passing an anonymous function like \code{function(x) weight_transcripts(x, type = 'equal')}.
#' @param filter_fn Optional filter function to filter the transcripts before plotting.
#' See the \code{filter_fn} parameter in \code{\link{filter_coverage}} for more details.
#' Providing a filter fucntion here is the same as providing it in \code{\link{filter_coverage}}
#' and then passing the result to this function.
#' @param detailed logical, if \code{TRUE}, also plot the top 10 transcripts with the highest
#' read counts for each length bin.
#' @return a ggplot2 object of the coverage plot(s)
#' @examples
#' # Create a BAM file with minimap2_realign
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
#' # Plot the coverages directly from the BAM file
#' plot_coverage(file.path(outdir, 'realign2transcript.bam'))
#'
#' # Get the coverage information first
#' coverage <- get_coverage(file.path(outdir, 'realign2transcript.bam')) |>
#'   dplyr::filter(read_counts > 2) |> # Filter out transcripts with read counts < 3
#'   filter_coverage(filter_fn = convolution_filter) # Filter out transcripts with sharp drops / rises
#' # Plot the filtered coverages
#' plot_coverage(coverage, detailed = TRUE)
#' # filtering function can also be passed directly to plot_coverage
#' plot_coverage(file.path(outdir, 'realign2transcript.bam'), filter_fn = convolution_filter)
#' @md
#' @export
plot_coverage <- function(x, quantiles = c(0, 0.2375, 0.475, 0.7125, 0.95, 1),
    length_bins = c(0, 1, 2, 5, 10, Inf), weight_fn = weight_transcripts,
    filter_fn, detailed = FALSE) {

  if (is.character(x) || is(x, "GAlignments")) {
    transcript_info <- get_coverage(x)
  } else {
    transcript_info <- x
  }

  if (!missing(filter_fn)) {
    transcript_info <- filter_coverage(transcript_info, filter_fn)
  }

  if (missing(length_bins)) {
    message("Using quantiles to bin transcripts.")
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
        min_length = round(min(tr_length) / 1000, 2),
        max_length = round(max(tr_length) / 1000, 2),
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
    dplyr::mutate(weight = weight_fn(read_counts)) |>
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
  if (detailed) {
    plot_list <- list(main = p)
    for (i in unique(transcript_info$length_bin)) {
      p_i <- transcript_info |>
        dplyr::filter(length_bin == i) |>
        dplyr::arrange(desc(read_counts)) |>
        head(10) |>
        tidyr::pivot_longer(paste0("coverage_", 1:100), names_to = "x", values_to = "coverage") |>
        dplyr::mutate(x = as.numeric(gsub("coverage_", "", x))) |>
        ggplot2::ggplot(aes(x = x, y = coverage, color = transcript, linewidth = read_counts)) +
        ggplot2::geom_line() +
        ggplot2::xlab("% Position in transcript") +
        ggplot2::ylab("Weighted average coverage") +
        ggplot2::theme_minimal() +
        ggplot2::guides(color = ggplot2::guide_legend(title = "Transcript ID"))
      plot_list[[i]] <- p_i
    }
    return(cowplot::plot_grid(plotlist = plot_list))
  } else {
    return(p)
  }
}

#' @importFrom GenomicAlignments seqnames
#' @importFrom tibble as_tibble
#' @importFrom dplyr rename filter pull
#' @importFrom GenomicRanges granges coverage
get_coverage_RleList <- function(bam, isoform = NULL) {
  if (!is.null(isoform)) {
    bam <- bam[GenomicAlignments::seqnames(bam) %in% isoform]
  }

  transcript_names <- table(GenomicAlignments::seqnames(bam)) |>
    as.data.frame() |>
    dplyr::rename(transcript = Var1, read_counts = Freq) |>
    tibble::as_tibble() |>
    dplyr::filter(read_counts > 0) |>
    dplyr::pull(transcript)

  cover <- bam |>
    GenomicRanges::granges() |>
    GenomicRanges::coverage()
  cover <- cover[transcript_names]
  return(cover)
}

#' @importFrom tibble as_tibble
sample_coverage <- function(coverage, read_counts_df) {
  read_counts <- read_counts_df$Freq[match(names(coverage), read_counts_df$Var1)]
  cover <- coverage |>
    sapply(function(x) {
      x[round(seq(1, length(x), length.out = 100), 0)] |>
        as.integer()
    }) |>
    t() |>
    as.data.frame()
  colnames(cover) <- paste0("coverage_", 1:100)
  cover <- cover / read_counts
  cover <- cover |>
    tibble::as_tibble(rownames = "transcript")
  cover$tr_length <- sapply(coverage, length)
  cover$read_counts <- read_counts
  return(cover)
}

#' @importFrom GenomicFeatures threeUTRsByTranscript fiveUTRsByTranscript
remove_UTR_coverage <- function(coverage_rlel, annotation) {
  transcript_names <- names(coverage_rlel)

  message("Getting UTRs...")
  threeUTR.data <- GenomicFeatures::threeUTRsByTranscript(annotation, use.names = T)
  fiveUTR.data <- GenomicFeatures::fiveUTRsByTranscript(annotation, use.names = T)
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
    lapply(seq_along(coverage_rlel), function(x) {
      coverage_rlel[[x]][(fiveUTR.len[x] + 1):(length(coverage_rlel[[x]]) - threeUTR.len[x])]
    })
  names(cover) <- transcript_names
  return(cover)
}
