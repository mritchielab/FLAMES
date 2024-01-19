#' Isoform identification
#' @description Long-read isoform identification with FLAMES or bambu.
#' @param annotation Path to annotation file. If configured to use bambu, the annotation
#' must be provided as GTF file.
#' @param genome_fa The file path to genome fasta file.
#' @param genome_bam File path to BAM alignment file. Multiple files could be provided.
#' @param outdir The path to directory to store all output files.
#' @param config Parsed FLAMES configurations.
#' @return The updated annotation and the transcriptome assembly will be saved in the
#' output folder as \code{isoform_annotated.gff3} (GTF if bambu is selected) and
#' \code{transcript_assembly.fa} respectively.
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
#'     config <- jsonlite::fromJSON(system.file("extdata/SIRV_config_default.json", package = "FLAMES"))
#'     minimap2_align(
#'         config = config,
#'         fa_file = genome_fa,
#'         fq_in = fastq1,
#'         annot = annotation,
#'         outdir = outdir
#'     )
#' \dontrun{
#'     find_isoform(
#'         annotation = annotation, genome_fa = genome_fa,
#'         genome_bam = file.path(outdir, "align2genome.bam"),
#'         outdir = outdir, config = config
#'     )
#' }
#' }
find_isoform <- function(annotation, genome_fa, genome_bam, outdir, config) {
    # pipeline types: singe_cell, single_cell_multisample, bulk
    cat(format(Sys.time(), "%X %a %b %d %Y"), "find_isoform\n")
    if (config$pipeline_parameters$bambu_isoform_identification) {
        find_isoform_bambu(annotation, genome_fa, genome_bam, outdir, config)
    } else {
        find_isoform_flames(annotation, genome_fa, genome_bam, outdir, config)
    }
}

#' @importFrom bambu writeToGTF prepareAnnotations bambu
#' @importFrom withr with_package
#' @importFrom SummarizedExperiment assays rowRanges
find_isoform_bambu <- function(annotation, genome_fa, genome_bam, outdir, config) {
    bambuAnnotations <- bambu::prepareAnnotations(annotation)
    # Tmp fix: remove withr if bambu imports seqlengths properly
    # https://github.com/GoekeLab/bambu/issues/255
    # min.readCount seems to cause errors
    # https://github.com/GoekeLab/bambu/issues/364

    bambu_out <- withr::with_package("GenomeInfoDb", 
        bambu::bambu(
                reads = genome_bam, 
                annotations = bambuAnnotations, 
                genome = genome_fa, 
                quant = TRUE, 
                discovery = TRUE,
                lowMemory = TRUE,
                NDR = config$isoform_parameters$bambu_ndr,
                ncore = ifelse(is.vector(genome_bam), length(genome_bam), 1)
        ))

    bambu::writeToGTF(SummarizedExperiment::rowRanges(bambu_out), file.path(outdir, "isoform_annotated_unfiltered.gtf")) 
    if (is.null(config$isoform_parameters$bambu_trust_reference) || config$isoform_parameters$bambu_trust_reference) {
        bambu_out <- bambu_out[base::rowSums(SummarizedExperiment::assays(bambu_out)$counts)>=1,]
    } else {
        bambu_out <- bambu_out[base::rowSums(SummarizedExperiment::assays(bambu_out)$counts)>=config$isoform_parameters$min_sup_cn,]
    }

    isoform_gtf <- file.path(outdir, "isoform_annotated.gtf") # Bambu outputs GTF
    bambu::writeToGTF(SummarizedExperiment::rowRanges(bambu_out), isoform_gtf) 
    annotation_to_fasta(isoform_gtf, genome_fa, outdir)

    # isoform_objects <- list(transcript_dict = NULL, transcript_dict_i = parse_gff_tree(isoform_gtf)$transcript_dict)
    # isoform_objects
}

#' @importFrom reticulate import_from_path
#' @importFrom Rsamtools indexFa
#' @importFrom basilisk basiliskRun
find_isoform_flames <- function(annotation, genome_fa, genome_bam, outdir, config) {
    if (length(genome_bam) == 1) {
        ret <- basiliskRun(env = flames_env, fun = function(gff3, genome, iso, tss, fa, tran, ds, conf, raw) {
            python_path <- system.file("python", package = "FLAMES")
            find <- reticulate::import_from_path("find_isoform", python_path)
            ret <- find$find_isoform(gff3, genome, iso, tss, fa, tran, ds, conf, raw)
            ret
        },
        gff3 = annotation, genome = genome_bam, iso = file.path(outdir, "isoform_annotated.gff3"), tss = file.path(outdir, "tss_tes.bedgraph"), fa = genome_fa, tran = file.path(outdir, "transcript_assembly.fa"), ds = config$isoform_parameters$downsample_ratio, conf = config, raw = ifelse(config$isoform_parameters$generate_raw_isoform, file.path(outdir, "splice_raw.gff3"), FALSE)
        )
    } else {
        ret <- basiliskRun(env = flames_env, fun = function(gff3, genome, iso, tss, fa, tran, ds, conf, raw) {
            python_path <- system.file("python", package = "FLAMES")
            find <- reticulate::import_from_path("find_isoform", python_path)
            ret <- find$find_isoform_multisample(gff3, genome, iso, tss, fa, tran, ds, conf, raw)
            ret
        },
        gff3 = annotation, genome = genome_bam, iso = file.path(outdir, "isoform_annotated.gff3"), tss = file.path(outdir, "tss_tes.bedgraph"), fa = genome_fa, tran = file.path(outdir, "transcript_assembly.fa"), ds = config$isoform_parameters$downsample_ratio, conf = config, raw = ifelse(config$isoform_parameters$generate_raw_isoform, file.path(outdir, "splice_raw.gff3"), FALSE)
        )
    }
    # we then need to use Rsamtools to index transcript_fa
    Rsamtools::indexFa(file.path(outdir, "transcript_assembly.fa")) # index the output fa file
}

#' GTF/GFF to FASTA conversion
#' @description convert the transcript annotation to transcriptome assembly as FASTA file. The
#' genome annotation is first imported as TxDb object and then used to extract transcript sequence
#' from the genome assembly.
#' @param isoform_annotation Path to the annotation file (GTF/GFF3)
#' @param genome_fa The file path to genome fasta file.
#' @param outdir The path to directory to store the transcriptome as \code{transcript_assembly.fa}.
#' @param extract_fn (optional) Function to extract \code{GRangesList} from the genome TxDb object.
#' E.g. \code{function(txdb){GenomicFeatures::cdsBy(txdb, by="tx", use.names=TRUE)}}
#' @return Path to the outputted transcriptome assembly
#'
#' @importFrom Biostrings readDNAStringSet writeXStringSet
#' @importFrom GenomicFeatures extractTranscriptSeqs makeTxDbFromGFF
#' @importFrom Rsamtools indexFa
#'
#' @examples
#' fasta <- annotation_to_fasta(system.file("extdata/rps24.gtf.gz", package = "FLAMES"), system.file("extdata/rps24.fa.gz", package = "FLAMES"), tempdir())
#' cat(readChar(fasta, nchars = 1e3))
#'
#' @export
annotation_to_fasta <- function(isoform_annotation, genome_fa, outdir, extract_fn) {
  # check if all the transcript in the annotation is stranded
  annotation_d <- read.csv(isoform_annotation, sep = "\t", 
                    header = FALSE, stringsAsFactors = FALSE, 
                    comment.char = "#")
  strands  <- annotation_d[,7]
  if (any(strands == '.')) {
        strands[strands == '.'] <- '+'
        annotation_d[,7] <- strands
        modified_gtf <- paste0(tempfile(),'/tmp.gtf')
        dir.create(dirname(modified_gtf))
        write.table(annotation_d, modified_gtf, sep = "\t", 
                    row.names = FALSE, quote = FALSE, col.names = FALSE)
        isoform_annotation <- modified_gtf
    }
  rm(annotation_d, strands)

  out_file <- file.path(outdir, "transcript_assembly.fa")

  dna_string_set <- Biostrings::readDNAStringSet(genome_fa)
  names(dna_string_set) <- gsub(" .*$", "", names(dna_string_set))
  txdb <- GenomicFeatures::makeTxDbFromGFF(isoform_annotation)
  if (missing(extract_fn)) {
    tr_string_set <- GenomicFeatures::extractTranscriptSeqs(dna_string_set, txdb,
      use.names = TRUE)
  } else {
    extracted_grl<- extract_fn(txdb)
    tr_string_set <- GenomicFeatures::extractTranscriptSeqs(dna_string_set, extracted_grl)
    # additional arguments are allowed only when 'transcripts' is not a GRangesList object
  }

  if (length(names(tr_string_set)) > length(unique(names(tr_string_set)))) {
    cat("Duplicated transcript IDs present, removing ...")
    tr_string_set <- tr_string_set[unique(names(tr_string_set))]
  }

  Biostrings::writeXStringSet(tr_string_set, out_file)
  Rsamtools::indexFa(out_file)

  return(out_file)
}

#' Parse FLAMES' GFF output
#' @description Parse FLAMES' GFF ouputs into a Genomic Ranges List
#' @param file the GFF file to parse
#' @return A Genomic Ranges List
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
#'     config <- jsonlite::fromJSON(system.file("extdata/SIRV_config_default.json", package = "FLAMES"))
#'     minimap2_align(
#'         config = config,
#'         fa_file = genome_fa,
#'         fq_in = fastq1,
#'         annot = annotation,
#'         outdir = outdir
#'     )
#'     find_isoform(
#'         annotation = annotation, genome_fa = genome_fa,
#'         genome_bam = file.path(outdir, "align2genome.bam"),
#'         outdir = outdir, config = config
#'     )
#'     grlist <- get_GRangesList(file = file.path(outdir, "isoform_annotated.gff3"))
#' }
#' @importFrom rtracklayer import
#' @export
get_GRangesList <- function(file) {
    isoform_gr <- rtracklayer::import(file, feature.type = c("exon", "utr"))
    if (grepl("\\.gff3(\\.gz)?$", file)) {
        isoform_gr$Parent <- as.character(isoform_gr$Parent)
        isoform_gr$transcript_id <- unlist(lapply(strsplit(isoform_gr$Parent, split = ":"), function(x) {
            x[2]
        }))
    }
    #    if (!is.null("gene")) {
    #        isoform_gr <- isoform_gr[isoform_gr$gene_id == gene]
    #    }
    isoform_grl <- S4Vectors::split(isoform_gr, isoform_gr$transcript_id)
    return(isoform_grl)
}
