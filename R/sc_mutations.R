variant_count_tb <- function(bam_path, seqname, pos, indel, barcodes) {
  # allele by barcode matrix (value: read count)
  tryCatch(
    {
      variant_count_matrix( # throws Rcpp::exception when no reads at pos
        bam_path = bam_path,
        seqname = seqname, pos = pos, indel = indel, barcodes = barcodes
      ) |>
        tibble::as_tibble(rownames = "allele") |>
        # pivot to long format: allele, barcode, allele_count
        tidyr::pivot_longer(
          cols = -tidyselect::matches("allele"),
          values_to = "allele_count", names_to = "barcode"
        ) |>
        dplyr::group_by(barcode) |>
        dplyr::mutate(
          cell_total_reads = sum(allele_count),
          pct = allele_count / cell_total_reads,
          pos = pos, seqname = seqname
        ) |>
        dplyr::ungroup()
    },
    error = function(e) {
      if (inherits(e, "Rcpp::exception") & conditionMessage(e) == "Failed to fetch an alignment") {
        message(paste0("No reads found at ", seqname, ":", pos, " in ", bam_path))
        message("Returning empty tibble")
        return(tibble::tibble())
      } else {
        stop(e)
      }
    }
  )
}

#' Variant count for single-cell data
#'
#' Count the number of reads supporting each variants at the given positions for each cell.
#'
#' @param bam_path character(1) or character(n): path to the bam file(s) aligned to the
#' reference genome (NOT the transcriptome! Unless the postions are also from the transcriptome).
#' @param seqnames character(n): chromosome names of the postions to count alleles.
#' @param positions integer(n): positions, 1-based, same length as seqnames. The positions to count alleles.
#' @param indel logical(1): whether to count indels (TRUE) or SNPs (FALSE).
#' @param barcodes character(n) when bam_path is a single file, or list of character(n)
#' when bam_path is a list of files paths. The cell barcodes to count alleles for.
#' Only reads with these barcodes will be counted.
#' @param threads integer(1): number of threads to use. Maximum number of threads is
#' the number of bam files * number of positions.
#' @return A tibble with columns: allele, barcode, allele_count, cell_total_reads, pct, pos, seqname.
#' @examples
#' outdir <- tempfile()
#' dir.create(outdir)
#' genome_fa <- file.path(outdir, "rps24.fa")
#' R.utils::gunzip(filename = system.file("extdata/rps24.fa.gz", package = "FLAMES"), destname = genome_fa, remove = FALSE)
#' download.file("https://raw.githubusercontent.com/mritchielab/FLAMES/devel/tests/testthat/demultiplexed.fq",
#'   destfile = file.path(outdir, "demultipelxed.fq")
#' ) # can't be bothered to run demultiplexing again
#' if (is.character(locate_minimap2_dir())) {
#'   minimap2_align( # align to genome
#'     config = jsonlite::fromJSON(system.file("extdata/SIRV_config_default.json", package = "FLAMES")),
#'     fa_file = genome_fa,
#'     fq_in = file.path(outdir, "demultipelxed.fq"),
#'     annot = system.file("extdata/rps24.gtf.gz", package = "FLAMES"),
#'     outdir = outdir
#'   )
#'   snps_tb <- sc_mutations(
#'     bam_path = file.path(outdir, "align2genome.bam"),
#'     seqnames = c("chr14", "chr14"),
#'     positions = c(1260, 2714), # positions of interest
#'     indel = FALSE,
#'     barcodes = read.delim(system.file("extdata/bc_allow.tsv.gz", package = "FLAMES"), header = FALSE)$V1
#'   )
#'   head(snps_tb)
#'   snps_tb |>
#'     dplyr::filter(pos == 1260) |>
#'     dplyr::group_by(allele) |>
#'     dplyr::summarise(count = sum(allele_count)) # should be identical to samtools pileup
#' }
#' @export
sc_mutations <- function(bam_path, seqnames, positions, indel = FALSE, barcodes, threads = 1) {
  stopifnot(
    "seqnames not the same length as positions" =
      length(seqnames) == length(positions)
  )

  if (length(bam_path) == 1) {
    # single bam file, parallelize over positions
    stopifnot(
      "barcodes must be a character vector" =
        is.character(barcodes)
    )
    message(paste0(format(Sys.time(), "%H:%M:%S "), "Got 1 bam file, parallelizing over each position ..."))
    variants <- parallel::mcmapply(
      FUN = function(seqname, pos) {
        variant_count_tb(bam_path, seqname, pos, indel, barcodes)
      },
      seqname = seqnames, pos = positions, SIMPLIFY = FALSE, mc.cores = threads
    )
  } else {
    # multiple bam files, parallelize over bam files
    stopifnot(
      "barcodes must be a list of character vectors, same length as bam_path" =
        is.list(barcodes) & length(barcodes) == length(bam_path)
    )
    # data frame of all combinations between (seqname, pos) and (bam_path, barcodes)
    args_grid <- expand.grid(
      mutation_index = seq_along(positions),
      bam_index = seq_along(bam_path),
      stringsAsFactors = FALSE
    ) |>
      dplyr::mutate(
        seqname = seqnames[mutation_index],
        pos = positions[mutation_index],
        sample_bam = bam_path[bam_index],
        sample_barcodes = barcodes[bam_index]
      ) |>
      dplyr::select(-mutation_index, -bam_index)

    message(paste0(format(Sys.time(), "%H:%M:%S "), "Multi-threading over bam files x positions ..."))
    variants <- parallel::mcmapply(
      FUN = function(sample_bam, seqname, pos, sample_barcodes) {
        variant_count_tb(sample_bam, seqname, pos, indel, sample_barcodes) |>
          dplyr::mutate(bam_file = sample_bam)
      },
      sample_bam = args_grid$sample_bam, seqname = args_grid$seqname,
      pos = args_grid$pos, sample_barcodes = args_grid$sample_barcodes,
      SIMPLIFY = FALSE, mc.cores = threads
    )
  }

  message(paste0(format(Sys.time(), "%H:%M:%S "), "Merging results ..."))
  variants <- variants |>
    dplyr::bind_rows()
  return(variants)
}


extract_nt <- function(ref, seqname, pos) {
  mapply(function(seqname, pos) {
    as.character(ref[[seqname]][pos])
  }, seqname, pos)
}

homopolymer_pct <- function(ref, seqname, pos, include_alt = FALSE, n = 3, threads = 1) {
  parallel::mcmapply(
    FUN = function(seqname, pos, include_alt) {
      if (pos == 1 | pos == length(ref[[seqname]])) {
        return(TRUE) # variant at the ends should not be considered
      }
      start <- max(1, pos - n)
      end <- min(length(ref[[seqname]]), pos + n)
      if (include_alt) {
        ref[[seqname]][start:end] |>
          as.character() |>
          strsplit("") |>
          unlist() |>
          table() |>
          as.data.frame() |>
          dplyr::mutate(pct = Freq / sum(Freq)) |>
          dplyr::pull(pct) |>
          max()
      } else {
        # exclude the position itself
        ref[[seqname]][c(start:(pos - 1), (pos + 1):end)] |>
          as.character() |>
          strsplit("") |>
          unlist() |>
          table() |>
          as.data.frame() |>
          dplyr::mutate(pct = Freq / sum(Freq)) |>
          dplyr::pull(pct) |>
          max()
      }
    }, seqname, pos, include_alt, mc.cores = threads
  )
}

# WIP: too much sequencing errrors / splice sites
# find variants in a single grange
find_variants_grange <- function(bam_path, reference, gene_grange, min_nucleotide_depth,
                                 names_from) {
  # read bam file
  mutations <- Rsamtools::pileup(bam_path,
    pileupParam = Rsamtools::PileupParam(
      max_depth = .Machine$integer.max - 1, min_base_quality = 0, min_mapq = 0,
      min_nucleotide_depth = min_nucleotide_depth, min_minor_allele_depth = 0,
      distinguish_strands = FALSE, distinguish_nucleotides = TRUE,
      ignore_query_Ns = TRUE, include_deletions = TRUE, include_insertions = TRUE,
      left_bins = NULL, query_bins = NULL, cycle_bins = NULL
    ),
    scanBamParam = Rsamtools::ScanBamParam(which = gene_grange)
  ) |>
    dplyr::select(-which_label) |>
    # move insertion postions to the base at which insertion start
    # like in IGV
    # don't count insertion in sum
    dplyr::mutate(
      pos = ifelse(nucleotide == "+", pos + 1, pos),
      counts_no_ins = ifelse(nucleotide != "+", count, 0)
    ) |>
    dplyr::group_by(seqnames, pos) |>
    dplyr::mutate(sum = sum(counts_no_ins)) |>
    dplyr::ungroup() |>
    dplyr::select(-counts_no_ins) |>
    dplyr::mutate(
      freq = count / sum,
      ref = factor(extract_nt(reference, seqnames, pos))
    ) |>
    dplyr::filter(as.character(nucleotide) != as.character(ref))

  if (nrow(mutations) == 0) {
    return(mutations)
  } else {
    mutations$bam_path <- bam_path
    if (names_from %in% colnames(S4Vectors::mcols(gene_grange))) {
      mutations$region <- S4Vectors::mcols(gene_grange)[, names_from]
    } else {
      mutations$region <- NA # no gene name / gap
    }
    return(mutations)
  }
}

#' bulk variant identification
#'
#' Treat each bam file as a bulk sample and identify variants against the reference
#'
#' Each bam file is treated as a bulk sample to perform pileup and identify variants.
#' You can run \code{sc_mutations} with the variants identified with this function
#' to get single-cell allele counts. Note that reference genome FASTA files may have
#' the chromosome names field as `>chr1 1` instead of `>chr1`. You may need to remove
#' the trailing number to match the chromosome names in the bam file, for example with
#' \code{names(ref) <- sapply(names(ref), function(x) strsplit(x, " ")[[1]][1])}.
#'
#' @param bam_path character(1) or character(n): path to the bam file(s) aligned to the
#' reference genome (NOT the transcriptome!).
#' @param reference DNAStringSet: the reference genome
#' @param annotation GRanges: the annotation of the reference genome. You can load
#' a GTF/GFF annotation file with \code{anno <- rtracklayer::import(file)}.
#' @param min_nucleotide_depth integer(1): minimum read depth for a position to be
#' considered a variant.
#' @param threads integer(1): number of threads to use. Threading is done over each
#' annotated region and (if \code{annotated_region_only = FALSE}) unannotated gaps for
#' each bam file.
#' @param homopolymer_window integer(1): the window size to calculate the homopolymer
#' percentage. The homopolymer percentage is calculated as the percentage of the most
#' frequent nucleotide in a window of \code{-homopolymer_window} to \code{homopolymer_window}
#' nucleotides around the variant position, excluding the variant position itself.
#' Calculation of the homopolymer percentage is skipped when \code{homopolymer_window = 0}.
#' This is useful for filtering out Nanopore sequencing errors in homopolymer regions.
#' @param annotated_region_only logical(1): whether to only consider variants outside
#' annotated regions. If \code{TRUE}, only variants outside annotated regions will be
#' returned. If \code{FALSE}, all variants will be returned, which could take significantly
#' longer time.
#' @param names_from character(1): the column name in the metadata column of the annotation
#' (\code(mcols(annotation)[, names_from])) to use for the \code{region} column in the output.
#' @return A tibble with columns: seqnames, pos, nucleotide, count, sum, freq, ref, region,
#' homopolymer_pct, bam_path The homopolymer percentage is calculated as the percentage of the
#' most frequent nucleotide in a window of \code{homopolymer_window} nucleotides around
#' the variant position, excluding the variant position itself.
#' @examples
#' outdir <- tempfile()
#' dir.create(outdir)
#' genome_fa <- file.path(outdir, "rps24.fa")
#' R.utils::gunzip(filename = system.file("extdata/rps24.fa.gz", package = "FLAMES"), destname = genome_fa, remove = FALSE)
#' download.file("https://raw.githubusercontent.com/mritchielab/FLAMES/devel/tests/testthat/demultiplexed.fq",
#'   destfile = file.path(outdir, "demultipelxed.fq")
#' ) # cant be bothered to run demultiplexing again
#' if (is.character(locate_minimap2_dir())) {
#'   minimap2_align( # align to genome
#'     config = jsonlite::fromJSON(system.file("extdata/SIRV_config_default.json", package = "FLAMES")),
#'     fa_file = genome_fa,
#'     fq_in = file.path(outdir, "demultipelxed.fq"),
#'     annot = system.file("extdata/rps24.gtf.gz", package = "FLAMES"),
#'     outdir = outdir
#'   )
#'   variants <- find_variants(
#'     bam_path = file.path(outdir, "align2genome.bam"),
#'     reference = genome_fa,
#'     annotation = GenomicRanges::GRanges("chr14", IRanges::IRanges(1, 1)),
#'     min_nucleotide_depth = 10
#'   )
#'   head(variants)
#' }
#' @export
find_variants <- function(bam_path, reference, annotation, min_nucleotide_depth = 100,
                          homopolymer_window = 3, annotated_region_only = FALSE,
                          names_from = "gene_name", threads = 1) {
  if (is.character(reference)) {
    message(paste0(format(Sys.time(), "%H:%M:%S "), "Reading reference ..."))
    reference <- Biostrings::readDNAStringSet(reference)
    # get rid of the `1` from >chr1 1
    names(reference) <- sapply(names(reference), function(x) strsplit(x, " ")[[1]][1])
  }
  if (is.character(annotation)) {
    message(paste0(format(Sys.time(), "%H:%M:%S "), "Reading annotation ..."))
    annotation <- rtracklayer::import(annotation) |>
      (\(x) x[S4Vectors::mcols(x)$type == "gene", ])()
  }

  if (!annotated_region_only) {
    message(paste0(format(Sys.time(), "%H:%M:%S "), "Adding unannotated gaps ..."))
    # GenomicRanges::gaps do not return the gap at the end of each chromosome,
    # so we add a dummy end to each sequence
    ends <- sapply(
      seq_along(reference),
      function(x) {
        GenomicRanges::GRanges(
          names(reference)[x],
          IRanges::IRanges(length(reference[[x]]), length(reference[[x]]))
        )
      }
    )
    # throws warning about not sequence levels in common, expected
    annotation <- c(annotation, do.call(c, ends))
    annotation <- c(annotation, GenomicRanges::gaps(annotation))
  }

  if (length(bam_path) == 1) {
    message(paste0(format(Sys.time(), "%H:%M:%S "), "Got 1 bam file, parallelizing over each region ..."))
    variants <- parallel::mclapply(
      sapply(seq_along(annotation), function(x) annotation[x]), function(grange) {
        find_variants_grange(bam_path, reference, grange, min_nucleotide_depth, names_from)
      },
      mc.cores = threads
    )
  } else {
    # multi-threading over bam files x granges
    message(paste0(format(Sys.time(), "%H:%M:%S "), "Got multiple bam files, preparing for multi-threading ..."))
    args_grid <- expand.grid(
      grange = seq_along(annotation), # does not work on GRanges directly
      bam = bam_path,
      stringsAsFactors = FALSE
    )

    message(paste0(format(Sys.time(), "%H:%M:%S "), "Multi-threading over bam files x ranges ..."))
    variants <- parallel::mcmapply(
      FUN = function(bam, grange) {
        find_variants_grange(bam, reference, annotation[grange], min_nucleotide_depth, names_from)
      },
      bam = args_grid$bam, grange = args_grid$grange, mc.cores = threads,
      SIMPLIFY = FALSE
    )
  }

  message(paste0(format(Sys.time(), "%H:%M:%S "), "Merging results ..."))
  variants <- dplyr::bind_rows(variants)
  if (nrow(variants) == 0) {
    # otherwise homopolymer_pct will fail
    return(variants)
  }

  if (homopolymer_window > 1) {
    message(paste0(format(Sys.time(), "%H:%M:%S "), "Calculating homopolymer percentages ..."))
    variants$homopolymer_pct <- homopolymer_pct(
      reference, variants$seqnames, variants$pos,
      include_alt = FALSE, n = homopolymer_window, threads = threads
    )
  }

  return(variants)
}
