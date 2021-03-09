#' Parse GFF3 file formats
#'
#' Parses GFF3 file formats into a usable list of four elements. Accepts \code{.gtf}, \code{.gff} formats as 
#' well as \code{.gz} compressed versions.
#' 
#' @param gff_file File path in gff, gtf.gz or gtf format to parse.
#'
#' @importFrom reticulate import_from_path
#' 
#' @examples 
#' gff3_parse <- parse_gff_tree(system.file("extdata/SIRV_anno.gtf", package="FLAMES"))
parse_gff_tree <- function(gff_file) {
    ret <- callBasilisk(flames_env, function(args) {
        python_path <- system.file("python", package="FlamesR")
        parse <- reticulate::import_from_path("parse_gene_anno", python_path)

        ret <- parse$parse_gff_tree(gff_file)
        names(ret) <- c("chr_to_gene", "transcript_dict", "gene_to_transcript", "transcript_to_exon")

        ret
    }, args=args)

    ret
}

