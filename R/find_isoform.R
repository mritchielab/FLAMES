#' @export
find_isoform <- function(annotation, genome_fa, genome_bam, outdir, config) {
    # pipeline types: singe_cell, single_cell_multisample, bulk
    if (config$pipeline_parameters$bambu_isoform_identification) {
        find_isoform_bambu(annotation, genome_fa, genome_bam, outdir, config)
    } else {
        find_isoform_flames(annotation, genome_fa, genome_bam, outdir, config)
    }
}

#' @importFrom bambu writeToGTF prepareAnnotations bambu
#' @importFrom withr with_package
find_isoform_bambu <- function(annotation, genome_fa, genome_bam, outdir, config) {
    bambuAnnotations <- bambu::prepareAnnotations(annotation)
    # Tmp fix: remove withr if bambu imports seqlengths properly
    # https://github.com/GoekeLab/bambu/issues/255
    bambu_out <- withr::with_package("GenomeInfoDb", bambu::bambu(reads = genome_bam, annotations = bambuAnnotations, genome = genome_fa, quant = FALSE))

    isoform_gtf <- file.path(outdir, "isoform_annotated.gtf") # Bambu outputs GTF
    bambu::writeToGTF(bambu_out, isoform_gtf) # bambu_out is the extended annotation
    annotation_to_fasta(isoform_gtf, genome_fa, outdir)

    # isoform_objects <- list(transcript_dict = NULL, transcript_dict_i = parse_gff_tree(isoform_gtf)$transcript_dict)
    # isoform_objects
}

#' @importFrom reticulate import_from_path
#' @importFrom Rsamtools indexFa
find_isoform_flames <- function(annotation, genome_fa, genome_bam, outdir, config) {
    if (length(genome_bam) == 1) {
        ret <- callBasilisk(flames_env, function(gff3, genome, iso, tss, fa, tran, ds, conf, raw) {
            python_path <- system.file("python", package = "FLAMES")
            find <- reticulate::import_from_path("find_isoform", python_path)
            ret <- find$find_isoform(gff3, genome, iso, tss, fa, tran, ds, conf, raw)
            ret
        },
        gff3 = annotation, genome = genome_bam, iso = file.path(outdir, "isoform_annotated.gff3"), tss = file.path(outdir, "tss_tes.bedgraph"), fa = genome_fa, tran = file.path(outdir, "transcript_assembly.fa"), ds = config$isoform_parameters$downsample_ratio, conf = config, raw = ifelse(config$isoform_parameters$generate_raw_isoform, file.path(outdir, "splice_raw.gff3"), FALSE)
        )
    } else {
        ret <- callBasilisk(flames_env, function(gff3, genome, iso, tss, fa, tran, ds, conf, raw) {
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

#' @importFrom Biostrings readDNAStringSet writeXStringSet
#' @importFrom GenomicFeatures extractTranscriptSeqs
#' @importFrom Rsamtools indexFa
#' @export
annotation_to_fasta <- function(isoform_annotation, genome_fa, outdir, out_file) {
    if (!missing(outdir) && !missing(out_file) && out_file != file.path(outdir, "transcript_assembly.fa")) {
        stop("Please specify only one of 'outdir' and 'out_file'.")
    }
    if (missing(out_file)) {
        out_file <- file.path(outdir, "transcript_assembly.fa")
    }

    dna_string_set <- Biostrings::readDNAStringSet(genome_fa)
    names(dna_string_set) <- gsub(" .*$", "", names(dna_string_set))
    tr_string_set <- GenomicFeatures::extractTranscriptSeqs(dna_string_set, get_GRangesList(isoform_annotation))
    Biostrings::writeXStringSet(tr_string_set, out_file)
    Rsamtools::indexFa(out_file)

    return(out_file)
}

get_GRangesList <- function(file) {
    isoform_gr <- rtracklayer::import(file, feature.type = c("exon", "utr"))
    if (grepl("\\.gff3$", file)) {
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
