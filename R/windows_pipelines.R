#' FLAMES Windows Bulk Pipeline
#' 
#' An implementation of the FLAMES pipeline designed to run on Windows, or any OS 
#' without access to minimap2, for read realignment.
#' This pipeline requires external read alignment, in betwen pipeline calls.
#' 
#' @details
#' This function, \code{bulk_windows_pipeline_setup} is the first step in the 3 step Windows FLAMES
#' bulk pipeline, and should be run first, read alignment undertaken,
#' then \code{windows_pipline_isoforms} should be run, read realignment performed, and finally
#' \code{windows_pipeline_quantification} should be run.
#' For each function, besides \code{bulk_windows_pipeline_setup}, a list \code{pipeline_variables} is returned, which
#' contains the information required to continue the pipeline. This list should be passed into each function, and updated
#' with the returned list. In the case of \code{bulk_windows_pipeline_setup}, \code{pipeline_variables} is
#' the list returned. See the vignette 'Vignette for FLAMES bulk on Windows' for more details.
#'
#' @param annot gene annotations file in gff3  format
#' @param fastq file path to input fastq file
#' @param outdir directory to store all output files.
#' @param genome_fa genome fasta file.
#' @param downsample_ratio downsampling ratio if performing downsampling analysis.
#' @param config_file JSON configuration file. If specified,  \code{config_file} overrides
#' all configuration parameters
#' 
#' @return a list \code{pipeline_variables} with the required variables for execution of later Windows pipeline
#' steps. File paths required to perform minimap2 alignment are given in pipeline_variables$return_files.
#' This list should be given as input for \code{windows_pipeline_isoforms} after minimap2 alignment has taken place; \code{windows_pipeline_isoforms} is the
#' continuation of this pipeline.
#' 
#' @example inst/examples/windows_bulk_pipeline.R
#' @export 
bulk_windows_pipeline_setup <- function(annot, fastq, in_bam=NULL, outdir, genome_fa,
                downsample_ratio=1, config_file) {
    # filenames for internal steps
    infq <- paste(outdir, "merged.fastq.gz", sep="/")
    bc_file <- paste(outdir, "pseudo_barcode_annotation.csv", sep="/")

    # create output directory if one doesn't exist
    if (!dir.exists(outdir)) {
        cat("Output directory does not exists: one is being created\n")
        dir.create(outdir)
        print(outdir)
    }

    if (is.null(in_bam)) {
        # use existing merge fastq if already exists
        if (file.exists(infq)) {
            cat(infq, " already exists, no need to merge fastq files\n")
        } else {
            # this preprocessing needs only be done if we are using a fastq_dir, instead
            # of a bam file for reads,
            cat("Preprocessing bulk fastqs...\n")
            # run the merge_bulk_fastq function as preprocessing
            merge_bulk_fastq(fastq, bc_file, infq)
        }
    } else {
        bc_file = NULL;
        fastq=NULL;
    }

    windows_pipeline_setup(annot, infq, in_bam, outdir, genome_fa, downsample_ratio, config_file, bc_file, bulk=TRUE)
}

#' Windows Single Cell FLAMES Pipeline
#' 
#' An implementation of the FLAMES pipeline designed to run on Windows, or any OS 
#' without access to minimap2, for read realignment.
#' This pipeline requires external read alignment, in betwen pipeline calls.
#' 
#' @details
#' This function, \code{sc_windows_pipeline_setup} is the first step in the 3 step Windows FLAMES
#' single cell pipeline, and should be run first, read alignment undertaken,
#' then \code{windows_pipline_isoforms} should be run, read realignment performed, and finally
#' \code{windows_pipeline_quantification} should be run.
#' For each function, besides \code{sc_windows_pipeline_setup}, a list \code{pipeline_variables} is returned, which
#' contains the information required to continue the pipeline. This list should be passed into each function, and updated
#' with the returned list. In the case of \code{sc_windows_pipeline_setup}, \code{pipeline_variables} is
#' the list returned. See the vignette 'Vignette for FLAMES bulk on Windows' for more details.
#'
#' @param annot gene annotations file in gff3  format
#' @param fastq file path to input fastq file
#' @param outdir directory to store all output files.
#' @param genome_fa genome fasta file.
#' @param downsample_ratio downsampling ratio if performing downsampling analysis.
#' @param config_file JSON configuration file. If specified,  \code{config_file} overrides
#' all configuration parameters
#' 
#' @return a list \code{pipeline_variables} with the required variables for execution of later Windows pipeline
#' steps. File paths required to perform minimap2 alignment are given in pipeline_variables$return_files.
#' This list should be given as input for \code{windows_pipeline_isoforms} after minimap2 alignment has taken place; \code{windows_pipeline_isoforms} is the
#' continuation of this pipeline.
#' @export 
sc_windows_pipeline_setup <- function(annot, fastq, in_bam=NULL, outdir, genome_fa,
                downsample_ratio=1, config_file) {
    if (is.null(in_bam)) {
        if (match_barcode) {
                if (!file.exists(reference_csv)) stop("reference_csv must exists.")
                infq <- paste(outdir, "matched_reads.fastq.gz", sep="/")
                bc_stat <- paste(outdir, "matched_barcode_stat", sep="/")
                match_cell_barcode(fastq, bc_stat, infq, reference_csv, MAX_DIST, UMI_LEN)
        } else infq = fastq
    }

    windows_pipeline_setup(annot, infq, in_bam, outdir, genome_fa, downsample_ratio, config_file, bulk=FALSE)
}
    

#' Hidden generic pipeline setup function
windows_pipeline_setup <- function(annot, fastq, in_bam=NULL, outdir, genome_fa,
                downsample_ratio=1, config_file,
                bc_file=NULL, bulk=FALSE) {

    if (!dir.exists(outdir)) {
        cat("Output directory does not exists: one is being created\n")
        dir.create(outdir)
        print(outdir)
    }
    if (is.null(config_file)) {
        stop("windows_pipeline requires a configuration file")
    } else {
        config = parse_json_config(config_file)
    }

    if (!config$pipeline_parameters$do_isoform_identification) {
        stop("Isoform Identification is required for FLAMES execution. Change this value in the configuration file or \
        set the argument as TRUE.")
    }

    cat("Running FLAMES pipeline...\n")
    # argument verificiation
    if (downsample_ratio > 1 || downsample_ratio <= 0) {
        stop("downsample_ratio should be between 0 and 1")
    }
    if (!file.exists(fastq)) stop(paste0("Make sure ", fastq, " exists."))
    if (!file.exists(annot)) stop(paste0("Make sure ", annot, " exists."))
    if (!file.exists(genome_fa)) stop(paste0("Make sure ", genome_fa, " exists."))

    using_bam = FALSE
    if (!is.null(in_bam)) {
        using_bam = TRUE
        bc_file = NULL
        fastq = NULL
        if (!file.exists(in_bam)) stop ("Make sure in_bam exists")
    } 

    pipeline_variables = list(
    # setup of internal arguments which hold output files and intermediate files
        isoform_gff3 = paste(outdir, "isoform_annotated.gff3", sep="/"),
        isoform_gff3_f = paste(outdir, "isoform_annotated.filtered.gff3", sep="/"),
        FSM_anno_out = paste(outdir, "isoform_FSM_annotation.csv", sep="/"),
        raw_splice_isoform = paste(outdir, "splice_raw.gff3", sep="/"),
        tss_tes_stat = paste(outdir, "tss_tes.bedgraph", sep="/"),
        transcript_fa = paste(outdir, "transcript_assembly.fa", sep="/"),
        transcript_fa_idx = paste(outdir, "transcript_assembly.fa.fai", sep="/"),
        tmp_bam = paste(outdir, "tmp_align.bam", sep="/"),
        tmp_bed = paste(outdir, "tmp_splice_anno.bed12", sep="/"),
        tmp_sam = paste(outdir, "tmp_align.sam", sep="/"),
        genome_bam = if (using_bam) in_bam else paste(outdir, "align2genome.bam", sep="/"),
        realign_bam = paste(outdir, "realign2transcript.bam", sep="/"),
        tr_cnt_csv = paste(outdir, "transcript_count.csv.gz", sep="/"),
        tr_badcov_cnt_csv = paste(outdir, "transcript_count.bad_coverage.csv.gz", sep="/"),
        annot=annot,
        fastq=fastq,
        in_bam=in_bam,
        outdir=outdir,
        genome_fa=genome_fa,
        downsample_ratio=downsample_ratio,
        config=config,
        bc_file=bc_file,
        using_bam=using_bam,
        return_files = NULL,
        is_bulk = bulk
    )

    cat("#### Input parameters:\n")
    print_config(config)
    cat("\tgene annotation:", annot, "\n")
    cat("\tgenome fasta:", genome_fa, "\n")
    if (using_bam) {
        cat("\tinput bam:", in_bam, "\n")
    } else cat("\tinput fastq:", fastq, "\n")
    cat("\toutput directory:", outdir, "\n")

    # align reads to genome
    if (!using_bam && config$pipeline_parameters$do_genome_alignment) {
        # at this point we need to stop the pipeline, return the files the user needs
        # to run minimap2, and then the pipeline can resume afterwards
        # to return for minimap2:
        # genome_fa, fastq, no_flank, bed12_junc
        pipeline_variables$return_files = list(
            genome_fa=genome_fa, 
            fastq=fastq, 
            no_flank=config$alignment_parameters$no_flank, 
            bed12_junc= if(config$alignment_parameters$use_junctions) pipeline_variables$tmp_bed else NULL,
            annot=annot,
            output_bam=pipeline_variables$genome_bam)
        return(pipeline_variables)
    } else {
        cat("#### Skip aligning reads to genome\n")
        windows_pipline_isoforms(pipline_variables)
    }
}

#' Windows Pipeline - Find Isoforms
#' 
#' This is the second step in the 3 step Windows FLAMES pipeline. 
#' Following this step, read realignment should be undertaken, using the file paths
#' given in the return pipeline_variables$return_files.
#' After this has been completed, the final pipeline step, \code{windows_pipeline_quantification} should be run,
#' giving the returned list from this function as input.
#' 
#' @param pipeline_variables the list returned from \code{windows_pipeline_isoforms}.
#' 
#' @return the updated \code{pipeline_variables} list, with information required for the final pipeline step.
#' 
#' @example inst/examples/windows_bulk_pipeline.R
#' @export
windows_pipeline_isoforms <- function(pipeline_variables) {
    # find isofrom
    cat("#### Read genne annotations\n")
    gff3_parse_result <- parse_gff_tree(pipeline_variables$annot)
    chr_to_gene = gff3_parse_result$chr_to_gene
    transcript_dict = gff3_parse_result$transcript_dict
    gene_to_transcript = gff3_parse_result$gene_to_transcript
    transcript_to_exon = gff3_parse_result$transcript_to_exon
    remove_similar_tr(gene_to_transcript, transcript_to_exon)

    cat("#### Find isoforms\n")
    transcript_to_junctions = list()
    for (tr in names(transcript_to_exon)) {
        transcript_to_junctions[[tr]] = blocks_to_junctions(transcript_to_exon[[tr]])
    }
    gene_dict <- get_gene_flat(gene_to_transcript, transcript_to_exon)
    chr_to_blocks <- get_gene_blocks(gene_dict, chr_to_gene, gene_to_transcript)
    group_bam2isoform(pipeline_variables$genome_bam, pipeline_variables$isoform_gff3, pipeline_variables$tss_tes_stat, 
            "", chr_to_blocks,
            gene_dict, transcript_to_junctions, 
            transcript_dict, pipeline_variables$genome_fa,
            config=pipeline_variables$config$isoform_parameters, 
            downsample_ratio=pipeline_variables$downsample_ratio,
            raw_gff3= 
            if (pipeline_variables$config$global_parameters$generate_raw_isoform) pipeline_variables$raw_splice_isoform else NULL
            )

    # get fasta
    isoform_gff3_parse <- parse_gff_tree(pipeline_variables$isoform_gff3)
    chr_to_gene_i <- isoform_gff3_parse$chr_to_gene
    transcript_dict_i <- isoform_gff3_parse$transcript_dict ## this is the only variable required after get_transcript_seq
    gene_to_transcript_i <- isoform_gff3_parse$gene_to_transcript
    transcript_to_exon_i <- isoform_gff3_parse$transcript_to_exon

    get_transcript_seq(pipeline_variables$genome_fa, pipeline_variables$transcript_fa, chr_to_gene_i, transcript_dict_i,
            gene_to_transcript_i, transcript_to_exon_i, ref_dict=if (pipeline_variables$config$realign_parameters$use_annotation) gff3_parse_result else NULL)

    pipeline_variables <- c(pipeline_variables, list(
        transcript_dict=transcript_dict,
        transcript_dict_i=transcript_dict_i
    ))
    # realign to transcript
    if (!pipeline_variables$using_bam && pipeline_variables$config$pipeline_parameters$do_read_realignment) {
        # the function needs to stop here again
        pipeline_variables$return_files <- list(
            transcript_fa=pipeline_variables$transcript_fa,
            fastq=pipeline_variables$fastq,
            out_bam=pipeline_variables$realign_bam
        )
        return (pipeline_variables)
    } else {
        cat("#### Skip read realignment\n")
        windows_pipeline_quantification(pipeline_variables)
    }
}

#' Windows Pipeline - Quantification
#' 
#' This is the final step in the 3 step Windows FLAMES pipeline. This should be run
#' after read realignment is performed, following \code{windows_pipeline_isoforms}.
#' 
#' @param pipeline_variables the list returned from \code{windows_pipeline_isoforms}, containing the information
#' required to perform the final step, quantification.
#' 
#' @return \code{windows_pipeline_quantification} returns a SummarizedExperiment object, or a SingleCellExperiment in the case
#' of this function being used for the FLAMES single cell pipeline, containing a count
#' matrix as an assay, gene annotations under metadata, as well as a list of the other
#' output files generated by the pipeline. The pipeline also outputs a number of output
#' files into the given \code{outdir} directory. These output files generated by the pipeline are:
#' \itemize{
#'  \item{transcript_count.csv.gz}{ - a transcript count matrix (also contained in the SummarizedExperiment)}
#'  \item{isoform_annotated.filtered.gff3}{ - isoforms in gff3 format (also contained in the SummarizedExperiment)}
#'  \item{transcript_assembly.fa}{ - transcript sequence from the isoforms}
#'  \item{align2genome.bam}{ - sorted BAM file with reads aligned to genome}
#'  \item{realign2transcript.bam}{ - sorted realigned BAM file using the transcript_assembly.fa as reference}
#'  \item{tss_tes.bedgraph}{ - TSS TES enrichment for all reads (for QC)}
#' }
#' 
#' @example inst/examples/bulk_windows_pipeline.R
#' @export
windows_pipeline_quantification <- function(pipeline_vars) {
    #quantification
    if (pipeline_vars$config$pipeline_parameters$do_transcript_quantification) {
        cat("#### Generating transcript count matrix\n")
        if (is.null(pipeline_vars$bc_file) || pipeline_vars$using_bam) {
            # sc_long_pipeline version of realigned bam
            parse_realign <- parse_realigned_bam(pipeline_vars$realign_bam, 
                pipeline_vars$transcript_fa_idx,
                pipeline_vars$config$isoform_parameters$Min_sup_cnt,
                pipeline_vars$config$transcript_counting$min_tr_coverage,
                pipeline_vars$config$transcript_counting$min_read_coverage)
        } else {
            parse_realign <- parse_realigned_bam(pipeline_vars$realign_bam, 
                pipeline_vars$transcript_fa_idx,
                pipeline_vars$config$isoform_parameters$Min_sup_cnt,
                pipeline_vars$config$transcript_counting$min_tr_coverage,
                pipeline_vars$config$transcript_counting$min_read_coverage,
                bc_file=pipeline_vars$bc_file)
        }
        tr_cnt = wrt_tr_to_csv(parse_realign$bc_tr_count_dict, 
            pipeline_vars$transcript_dict_i, pipeline_vars$tr_cnt_csv, 
            pipeline_vars$transcript_dict, 
            pipeline_vars$config$global_parameters$has_UMI)
        wrt_tr_to_csv(parse_realign$bc_tr_badcov_count_dict, 
            pipeline_vars$transcript_dict_i, 
            pipeline_vars$tr_badcov_cnt_csv, 
            pipeline_vars$transcript_dict, 
            pipeline_vars$config$global_parameters$has_UMI)
        annotate_filter_gff(pipeline_vars$isoform_gff3, 
            pipeline_vars$annot, 
            pipeline_vars$isoform_gff3_f, 
            pipeline_vars$FSM_anno_out, 
            tr_cnt, 
            pipeline_vars$config$isoform_parameters$Min_sup_cnt)
    } else {
        cat("#### Skip transcript quantification\n")
    }

    if (pipeline_vars$is_bulk) {
        return(generate_bulk_summarized(pipeline_vars$outdir))
    } else {
        return(generate_sc_singlecell(pipeline_vars$outdir))
    }
}
            

