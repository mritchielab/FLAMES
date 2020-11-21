#' Title
#'
#' DESC
#' 
#' @param gff_file File path in gff, gtf.gz or gtf format to parse.
#'
#' @importFrom reticulate import_from_path
#' @export
#'
#' @value List object with four elements: chr_to_gene, transcript_dict, gene_to_transcript, transcript_to_exon
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

