variant_count_tb <- function(bam_path, seqname, pos, indel, barcodes) {
  # allele by barcode matrix (value: read count)
  variant_count_matrix(
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
#' download.file('https://raw.githubusercontent.com/mritchielab/FLAMES/devel/tests/testthat/demultiplexed.fq', 
#'   destfile = file.path(outdir, 'demultipelxed.fq')) # can't be bothered to run demultiplexing again
#' if (is.character(locate_minimap2_dir())) {
#'   minimap2_align( # align to genome
#'       config = jsonlite::fromJSON(system.file('extdata/SIRV_config_default.json', package = 'FLAMES')),
#'       fa_file = genome_fa,
#'       fq_in = file.path(outdir, 'demultipelxed.fq'),
#'       annot = system.file("extdata/rps24.gtf.gz", package = "FLAMES"),
#'       outdir = outdir
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

  errors <- variants[sapply(variants, \(x) inherits(x, "try-error"))]
  if (length(errors) > 0) {
    warning(paste0(
      length(errors), " errors encountered out of ", length(variants), " positions * BAMs checked:"
    ))
    sapply(errors, \(x) x[1]) |>
      table() |>
      print()
  }

  variants <- variants[sapply(variants, \(x) !inherits(x, "try-error"))] |>
    dplyr::bind_rows()

  return(variants)
}


extract_nt <- function(ref, seqname, pos) {
  mapply(function(seqname, pos) {
    as.character(ref[[seqname]][pos])
  }, seqname, pos)
}

homopolymer_pct <- function(ref, seqname, pos, include_alt, n = 5) {
  mapply(function(seqname, pos, include_alt) {
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
   }, seqname, pos, include_alt)
}

# WIP: too much sequencing errrors / splice sites
find_variants <- function(bam_path, reference, gene_grange) {
  # read bam file
  pile <- Rsamtools::pileup(bam_path,
    pileupParam = Rsamtools::PileupParam(
      max_depth = .Machine$integer.max - 1, min_base_quality = 0, min_mapq = 0,
      min_nucleotide_depth = 5, min_minor_allele_depth = 0,
      distinguish_strands = FALSE, distinguish_nucleotides = TRUE,
      ignore_query_Ns = TRUE, include_deletions = TRUE, include_insertions = TRUE,
      left_bins = NULL, query_bins = NULL, cycle_bins = NULL
    ),
    scanBamParam = Rsamtools::ScanBamParam(which = gene_grange)
  )

  pile$ref <- extract_nt(reference, pile$seqnames, pile$pos) # add reference allele

  mutations <- pile |>
    dplyr::filter(as.character(nucleotide) != ref)

  # add variant frequency
  mutations$sum <- mapply(function(i, j) { # seqname, pos
    pile |>
      dplyr::filter(seqnames == i, pos == j) |>
      dplyr::pull(count) |>
      sum()
  }, mutations$seqnames, mutations$pos)
  mutations <- mutations |>
    # + counted twice for insertions (one for the current position and one for insertion)
    # still get freq > 1 (probably) due to ignore_query_Ns and min depths
    dplyr::mutate(freq = ifelse(nucleotide == "+" & count * 2 < sum, count * 2 / sum, count / sum))

  if (nrow(mutations) == 0) {
    return(mutations)
  } else {
    mutations$gene <- gene_grange$gene_name
    return(mutations)
  }
}
