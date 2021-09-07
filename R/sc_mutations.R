#'
#' FLAMES variant calling
#' 
#' Candidate SNVs identified with filtering by coverage threshold,
#' and allele frequency range. 
#' 
#' Takes the \code{SingleCellExperiment} object from \code{sc_long_pipeline} and the cell barcodes as \code{barcode}.
#' Alternatively, input can also be provided as \code{out_dir}, \code{genome_fa}, \code{annot}, \code{barcode}.
#' 
#' @param sce The \code{SingleCellExperiment} object from \code{sc_long_pipeline}
#' 
#' @param genome_fa (Optional) Reference genome FASTA file. Use this parameter is if you do not wish \code{sc_mutation} to use the
#' reference genome FASTA file from the \code{sce}'s metadata.
#' 
#' @param bam_short (Optional) short read alignment BAM file
#' 
#' @param out_dir (Optional) Output folder of sc_long_pipeline. Output files from this function will also be saved here.
#' Use this parameter if you do not have the \code{SingleCellExperiment} object.
#' 
#' @param barcode TSV file for cell barcodes
#' 
#' @param annot (Optional) The file path to gene annotation file in gff3 format. If provided as \code{FALSE} then the \code{isoform_annotated.gff3} from \code{sc_longread_pipeline} will be used, if not provided then the path in the \code{SingleCellExperiment} object will be used.
#' 
#' @param known_positions (Optional) A list of known positions, with by chromosome name followed by the position, e.g. ('chr1', 123, 'chr1', 124, 'chrX', 567)
#' 
#' @param min_cov The coverage threshod for filtering candidate SNVs
#' 
#' @param report_pct The allele frequency range for filtering candidate SNVs
#' 
#' @return
#' a \code{data.frame} containing the following columns:
#' \itemize{
#'  \item{chr}{ - the chromosome where the mutation is located}
#'  \item{position}
#'  \item{REF}{ - the reference allele}
#'  \item{ALT}{ - the alternative allele}
#'  \item{REF_frequency}{ - reference allele frequency}
#'  \item{REF_frequency_in_short_reads}{ - reference allele frequency in short reads (-1 when short reads not provided)}
#'  \item{hypergeom_test_p_value}
#'  \item{sequence_entrophy}
#'  \item{INDEL_frequency}
#'  \item{adj_p}{ - the adjusted p-value (by Benjamini–Hochberg correction)}
#' }
#' The table is sorted by decreasing adjusted P value.
#' 
#' files saved to out_dir/mutation:
#' \itemize{
#'      \item{ref_cnt.csv.gz}
#'      \item{alt_cnt.csv.gz}
#'      \item{allele_stat.csv.gz}
#'      \item{freq_summary.csv}
#' }
#' 
#' @importFrom reticulate import_from_path
#' @importFrom stats p.adjust
#' 
#' @export
sc_mutations <- function(sce, barcode_tsv, bam_short, out_dir, genome_fa, annot, known_positions=NULL, min_cov=100, report_pct=c(0.15, 0.85), test_mode=FALSE) {

    cat('\n\n\n\n\n\n')
    if (!missing(sce)) {
        if (missing(annot)) {
            annot <- sce@mdata$AnnotationFile
        }
        if (missing(genome_fa)) {
            genome_fa <- sce@mdata$genome_fa
        }
        if (missing(out_dir)) {
            out_dir <- sce@mdata$outdir
        }
    }

    python_path <- system.file("python", package = "FLAMES")
    callBasilisk(flames_env, function(fa, bam, dir, barcode, gff, positions, mincov, reportpct, test) {
        convert <- reticulate::import_from_path("sc_mutations", python_path)
        convert$sc_mutations(fa, bam, dir, barcode, gff, positions, mincov, reportpct, test)
    }, fa = genome_fa, bam = bam_short, dir = out_dir, barcode = barcode_tsv, gff=annot, positions = known_positions, mincov=min_cov, reportpct=report_pct, test = test_mode)

    allele_stat_csv <- file.path(out_dir, 'mutation', 'allele_stat.csv.gz')
    table <- read.csv(allele_stat_csv) 
    table$adj_p <- stats::p.adjust(table$hypergeom_test_p_value)
    table <- table[order(table$adj_p), ]
    rownames(table) <- NULL
    write.csv(table, file=gzfile(allele_stat_csv), row.names=FALSE)
    print.data.frame(head(table, n=20L))
    table
}
