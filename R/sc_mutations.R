#'
#' FLAMES variant calling
#' 
#' Candidate SNVs identified with filtering by coverage threshold,
#' and allele frequency range. 
#' 
#' @param fa_f Reference genome FASTA file
#' 
#' @param bam_short Optional short read alignment BAM file
#' 
#' @param out_dir Output folder of sc_long_pipeline. Output files from this function will also be saved here.
#' 
#' @param barcode TSV file for cell barcodes
#' 
#' @param gff_f Reference genome annotation file, in \code{GFF} or \code{GTF} format. If left with \code{NONE} then the \code{isoform_annotated.gff3} from \code{sc_longread_pipeline} will be used.
#' 
#' @param known_positions A list of known positions, with by chromosome name followed by the position, e.g. ('chr1', 123, 'chr1', 124, 'chrX', 567)
#' 
#' @param min_cov The coverage threshold for filtering candidate SNVs
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
sc_mutations <- function(fa_f, bam_short, out_dir, barcode_tsv, gff_f=NULL, known_positions=NULL, min_cov=100, report_pct=c(0.15, 0.85), test_mode=FALSE) {
    cat('\n\n\n\n\n\n')
    python_path <- system.file("python", package = "FLAMES")
    callBasilisk(flames_env, function(fa, bam, dir, barcode, gff, positions, mincov, reportpct, test) {
        convert <- reticulate::import_from_path("sc_mutations", python_path)
        convert$sc_mutations(fa, bam, dir, barcode, gff, positions, mincov, reportpct, test)
    }, fa = fa_f, bam = bam_short, dir = out_dir, barcode = barcode_tsv, gff=gff_f, positions = known_positions, mincov=min_cov, reportpct=report_pct, test = test_mode)

    allele_stat_csv <- file.path(out_dir, 'mutation', 'allele_stat.csv.gz')
    table <- read.csv(allele_stat_csv) 
    table$adj_p <- stats::p.adjust(table$hypergeom_test_p_value)
    table <- table[order(table$adj_p), ]
    rownames(table) <- NULL
    write.csv(table, file=gzfile(allele_stat_csv), row.names=FALSE)
    print.data.frame(head(table, n=20L))
    table
}
