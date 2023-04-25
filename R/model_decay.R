#' filter annotation for plotting coverages
#' 
#' @description Removes isoform annotations that could produce ambigious reads, such as isoforms
#' that only differ by the 5' / 3' end. This could be useful for plotting average
#' coverage plots.
#' 
#' @importFrom rtracklayer import
#' @importFrom S4Vectors split
#' @importFrom GenomicRanges strand
#' @importFrom BiocGenerics start end
#' @importFrom arrangements combinations
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
#' filtered_annotation <- filter_annotation(system.file('extdata/rps24.gtf.gz', package = 'FLAMES'), keep = 'tes_differ')
#' filtered_annotation
#'
#' @md
#' @export
filter_annotation <- function(annotation, keep = "tss_differ") {
  if (is.character(annotation)) {
    annotation <- rtracklayer::import(annotation, feature.type = "exon")
  }
  annotation <- S4Vectors::split(annotation, annotation$gene_id)

  min_differ <- function(x, y) {
    return(min(BiocGenerics::start(x)) != min(BiocGenerics::start(y)))
  }

  max_differ <- function(x, y) {
    return(max(BiocGenerics::end(x)) != max(BiocGenerics::end(y)))
  }

  tes_differ <- function(x, y) {
    if (as.logical(GenomicRanges::strand(x)[1] != GenomicRanges::strand(y)[1])) {
      return(TRUE)
    } else if (as.logical(GenomicRanges::strand(x)[1] == "+")) {
      return(max_differ(x, y))
    } else if (as.logical(GenomicRanges::strand(x)[1] == "-")) {
      return(min_differ(x, y))
    } else {
      return(min_differ(x, y) && max_differ(x, y))
    }
  }

  tss_differ <- function(x, y) {
    if (as.logical(GenomicRanges::strand(x)[1] != GenomicRanges::strand(y)[1])) {
      return(TRUE)
    } else if (as.logical(GenomicRanges::strand(x)[1] == "-")) {
      return(max_differ(x, y))
    } else if (as.logical(GenomicRanges::strand(x)[1] == "+")) {
      return(min_differ(x, y))
    } else {
      return(min_differ(x, y) && max_differ(x, y))
    }
  }

  check_fn <- switch(keep, single_transcripts = function(...) {
    stop("Line should not be reached")
  }, tss_differ = tss_differ, tes_differ = tes_differ, both = function(x, y) {
    tss_differ(x, y) && tes_differ(x, y)
  })

  filtered_annotation <- NULL
  for (gene in seq_along(annotation)) {
    transcripts <- S4Vectors::split(annotation[[gene]], annotation[[gene]]$transcript_id)
    if (length(transcripts) == 1) {
      filtered_annotation <- append(filtered_annotation, unlist(transcripts))
    } else if (keep == "single_transcripts") {
      next

      # remove isoforms that could spawn ambigious reads
    } else {
      pairwise_checks <- arrangements::combinations(x = names(transcripts),
        k = 2, replace = FALSE, layout = "row")
      kept_transcripts <- rep(TRUE, length(transcripts))
      names(kept_transcripts) <- names(transcripts)
      for (check in seq_len(nrow(pairwise_checks))) {
        i <- transcripts[[pairwise_checks[check, 1]]]
        j <- transcripts[[pairwise_checks[check, 2]]]
        if (!check_fn(i, j)) {
          kept_transcripts[pairwise_checks[check, ]] <- FALSE
        }
      }
      filtered_annotation <- append(filtered_annotation, unlist(transcripts[kept_transcripts]))
    }

  }
  return(filtered_annotation)
}

#' plot read coverages
#' 
#' @description Plot the average read coverages for each length bin or a 
#' perticular isoform
#' 
#' @importFrom GenomicFeatures makeTxDbFromGFF transcripts
#' @importFrom GenomicAlignments readGAlignments seqnames 
#' @importFrom GenomicRanges width strand granges coverage
#' @importFrom Rsamtools ScanBamParam
#' @importFrom tidyr as_tibble pivot_longer
#' @importFrom dplyr filter mutate group_by summarize_at summarise across
#' @importFrom ggplot2 ggplot geom_line aes
#' @importFrom stats weighted.mean
#'
#' @param annotation path to the GTF annotation file, or the parsed GenomicRanges
#' object.
#' @param isoform string vector, provide isoform names to plot the coverage for the
#' corresponding isoforms, or provide NULL to plot average coverages for each 
#' length bin
#' @param length_bins, numeric vector to specify the sizes to bin the isoforms by
#' @param bam, path to the BAM file (aligning reads to the transcriptome), or
#' the (GenomicAlignments::readGAlignments) parsed GAlignments object
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
#' if (is.character(locate_minimap2_dir())) {
#'     fasta <- annotation_to_fasta(annotation, genome_fa, outdir)
#'     minimap2_realign(
#'         config = jsonlite::fromJSON(system.file('extdata/SIRV_config_default.json', package = 'FLAMES')),
#'         fq_in = fastq1,
#'         outdir = outdir
#'     )
#'   plot_coverage(annotation = annotation, bam = file.path(outdir, 'realign2transcript.bam'))
#' }
#' @md
#' @export
plot_coverage <- function(annotation, bam, isoform = NULL, length_bins = c(0, 1,
  2, 5, 10, Inf)) {
  if (is.character(annotation)) {
    annotation <- annotation |>
      GenomicFeatures::makeTxDbFromGFF() |>
      GenomicFeatures::transcripts()
  } else {
    annotation <- annotation |>
      GenomicFeatures::makeTxDbFromGRanges() |>
      GenomicFeatures::transcripts()
  }

  if (!is(bam, "GAlignments")) {
    bam <- GenomicAlignments::readGAlignments(bam, param = Rsamtools::ScanBamParam(mapqFilter = 5))
  }

  read_counts <- table(GenomicAlignments::seqnames(bam))
  transcript_names <- intersect(annotation$tx_name, names(read_counts))
  annotation <- annotation[match(transcript_names, annotation$tx_name)]
  transcript_info <- data.frame(tr_length = GenomicRanges::width(annotation), read_counts = as.data.frame(read_counts[transcript_names])$Freq,
    strand = GenomicRanges::strand(annotation))
  transcript_info$length_bin <- cut(transcript_info$tr_length/1000, length_bins)

  cover <- bam |>
    GenomicRanges::granges() |>
    GenomicRanges::coverage() |>
    sapply(function(x) {
      x[round(seq(1, length(x), length.out = 100), 0)] |>
        as.integer()
    }) |>
    subset(select = transcript_names) |>
    t() |>
    as.data.frame()
  colnames(cover) <- paste0("coverage_", 1:100)
  cover <- cover[transcript_names, ]

  weight_covergae <- function(mat, read_counts) {
    sigmoid <- function(x) {
      exp(x)/(exp(x) + 1)
    }
    sigmoid((read_counts - mean(read_counts))/100)
  }

  cover <- cover/transcript_info$read_counts  # scale by read counts
  transcript_info <- cbind(transcript_info, cover)
  transcript_info$weight <- weight_covergae(mat, transcript_info$read_counts)
  if (!is.null(isoform)) {
    p <- transcript_info |>
      tidyr::as_tibble(rownames = "transcript") |>
      dplyr::filter(transcript %in% isoform) |>
      tidyr::pivot_longer(paste0("coverage_", 1:100), names_to = "x", values_to = "coverage") |>
      dplyr::mutate(x = as.numeric(gsub("coverage_", "", x))) |>
      ggplot2::ggplot(aes(x = x, y = coverage, color = transcript)) + geom_line()
    return(p)
  }

  mean_coverage <- transcript_info |>
    dplyr::group_by(length_bin) |>
    dplyr::summarise(dplyr::across(paste0("coverage_", 1:100), ~stats::weighted.mean(.,
      w = read_counts)))

  p <- mean_coverage |>
    tidyr::pivot_longer(paste0("coverage_", 1:100), names_to = "x", values_to = "coverage") |>
    dplyr::mutate(x = as.numeric(gsub("coverage_", "", x))) |>
    ggplot2::ggplot(aes(x = x, y = coverage, color = length_bin)) + geom_line()

  return(p)
}
