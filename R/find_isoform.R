#' @export
find_isoform <- function(annotation, genome_fa, genome_bam, outdir, config) {
    if (config$pipeline_parameters$bambu_isoform_identification) {
        isoform_objects <- find_isoform_bambu(annotation, genome_fa, genome_bam, outdir, config)
    } else {
        isoform_objects <- find_isoform_flames(annotation, genome_fa, genome_bam, outdir, config)
    }
    isoform_objects
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

    # Create transcriptome assembly .fa
    dna_string_set <- Biostrings::readDNAStringSet(genome_fa)
    names(dna_string_set) <- gsub(" .*$", "", names(dna_string_set))
    tr_string_set <- GenomicFeatures::extractTranscriptSeqs(dna_string_set, get_GRangesList(isoform_gtf))
    Biostrings::writeXStringSet(tr_string_set, file.path(outdir, "transcript_assembly.fa"))

    Rsamtools::indexFa(file.path(outdir, "transcript_assembly.fa"))
    # Todo: convert bambu_out (GRangesList) to transcript_dict directly
    isoform_objects <- list(transcript_dict = NULL, transcript_dict_i = parse_gff_tree(isoform_gtf)$transcript_dict)
    isoform_objects
}

#' @importFrom reticulate import_from_path
#' @importFrom Rsamtools indexFa
find_isoform_flames <- function(annotation, genome_fa, genome_bam, outdir, config) {
    ret <- callBasilisk(flames_env, function(gff3, genome, iso, tss, fa, tran, ds, conf, raw, seed) {
        python_path <- system.file("python", package = "FLAMES")
        find <- reticulate::import_from_path("find_isoform", python_path)
        ret <- find$find_isoform(gff3, genome, iso, tss, fa, tran, ds, conf, raw, seed)
        ret
    },
    gff3 = annotation, genome = genome_bam, iso = file.path(outdir, "isoform_annotated.gff3"), tss = file.path(outdir, "tss_tes.bedgraph"), fa = genome_fa, tran = file.path(outdir, "transcript_assembly.fa"), ds = config$isoform_parameters$downsample_ratio, conf = config, raw = ifelse(config$isoform_parameters$generate_raw_isoform, file.path(outdir, "splice_raw.gff3"), FALSE), seed = config$pipeline_parameters$seed
    )

    # we then need to use Rsamtools to index transcript_fa
    Rsamtools::indexFa(file.path(outdir, "transcript_assembly.fa")) # index the output fa file

    ret
}

find_isoform_flames_multisample <- function(gff3,
                                            genome_bams,
                                            isoform_gff3,
                                            tss_tes_stat,
                                            genome_fa,
                                            transcript_fa,
                                            downsample_ratio,
                                            config,
                                            raw) {
    ret <- callBasilisk(flames_env, function(gff3,
                                             genomes,
                                             iso,
                                             tss,
                                             fa,
                                             tran,
                                             ds,
                                             conf,
                                             raw) {
        python_path <- system.file("python", package = "FLAMES")

        find <-
            reticulate::import_from_path("find_isoform", python_path)
        ret <-
            find$find_isoform_multisample(gff3, genomes, iso, tss, fa, tran, ds, conf, raw)

        ret
    },
    gff3 = gff3, genomes = genome_bams, iso = isoform_gff3, tss = file.path(outdir, "tss_tes.bedgraph"),
    fa = genome_fa, tran = transcript_fa, ds = downsample_ratio, conf =
        config, raw = raw
    )

    # we then need to use Rsamtools to index transcript_fa
    Rsamtools::indexFa(transcript_fa) # index the output fa file

    ret
}
