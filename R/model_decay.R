#' filter annotation for plotting coverages
#' 
#' Removes isoform annotations that could produce ambigious reads, such as isoforms
#' that only differ by the 5' / 3' end. This could be useful for plotting average
#' coverage plots.
#' 
#' @param annotation path to the GTF annotation file, or the parsed GenomicRanges
#' object.
#' @param keep string, one of 'tss_differ' (only keep isoforms that all differ
#' by the transcription start site position), 'tes_differ' (only keep those that
#' differ by the transcription end site position), 'both' (only keep those that
#' differ by both the start and end site), or 'single_transcripts' (only keep
#' genes that contains a sinlge transcript).

filter_annotation <- function(annotation, keep = "tss_differ") {
  if (is.character(annotation)) {
    annotation <- rtracklayer::import(annotation, feature.type = "exon")
  }
  annotation <- S4Vectors::split(annotation, annotation$gene_id)

  min_differ <- function(x, y) {
    return(min(start(x)) != min(start(y)))
  }

  max_differ <- function(x, y) {
    return(max(end(x)) != max(end(y)))
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

plot_coverage <- function(annotation, isoform = NULL) {
  annotation <- annotation |>
    GenomicFeatures::makeTxDbFromGRanges() |>
    GenomicFeatures::transcripts()

  alig <- GenomicAlignments::readGAlignments("alignment/bp02_restranded.bam", param = Rsamtools::ScanBamParam(mapqFilter = 5))

  read_counts <- table(GenomicAlignments::seqnames(alig))
  transcript_names <- intersect(annotation$tx_name, names(read_counts))
  annotation <- annotation[match(transcript_names, annotation$tx_name)]
  transcript_info <- data.frame(tr_length = GenomicRanges::width(annotation), read_counts = as.data.frame(read_counts[transcript_names])$Freq,
    strand = GenomicRanges::strand(annotation))
  length_bins <- c(0, 1, 2, 5, 10, Inf)
  transcript_info$length_bin <- cut(transcript_info$tr_length/1000, length_bins)

  cover <- alig |>
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

  transcript_info <- cbind(transcript_info, cover[transcript_names, ])
  mean_coverage <- transcript_info |>
    dplyr::group_by(length_bin) |>
    dplyr::summarize_at(paste0("coverage_", 1:100), mean)

  mean_coverage |>
    tidyr::pivot_longer(paste0("coverage_", 1:100), names_to = "x", values_to = "coverage") |>
    dplyr::mutate(x = as.numeric(gsub("coverage_", "", x))) |>
    ggplot(aes(x = x, y = coverage, color = length_bin)) + geom_line()

  p <- transcript_info |>
    tidyr::as_tibble(rownames = "transcript") |>
    dplyr::filter(length_bin == "(1,2]") |>
    tidyr::pivot_longer(paste0("coverage_", 1:100), names_to = "x", values_to = "coverage") |>
    dplyr::mutate(x = as.numeric(gsub("coverage_", "", x))) |>
    ggplot(aes(x = x, y = coverage, color = transcript)) + geom_line()

  return(p)
}
