#' Parse Gff3 file
#'
#' @description
#' Parse a Gff3 file into 3 components: chromasome to gene name, a transcript dictionary, a gene to transcript dictionary
#' and a transcript to exon dictionary.
#' These components are returned in a named list.
#'
#' @param gff_file the file path to the gff3 file to parse
#' @return a named list with the elements
#' "chr_to_gene", "transcript_dict", "gene_to_transcript", "transcript_to_exon", containing
#' the data parsed from the gff3 file.
#' @importFrom reticulate import_from_path
#'
#' @examples
#' temp_path <- tempfile()
#' bfc <- BiocFileCache::BiocFileCache(temp_path, ask = FALSE)
#' file_url <-
#'     "https://raw.githubusercontent.com/OliverVoogd/FLAMESData/master/data"
#' gff <- bfc[[names(BiocFileCache::bfcadd(bfc, "GFF", paste(file_url, "SIRV_isoforms_multi-fasta-annotation_C_170612a.gtf", sep = "/")))]]
#' \donttest{
#' parsed_gff <- parse_gff_tree(gff)
#' }
#' @export
parse_gff_tree <- function(gff_file) {
    ret <- callBasilisk(flames_env, function(args) {
        python_path <- system.file("python", package = "FLAMES")
        parse <-
            reticulate::import_from_path("parse_gene_anno", python_path)

        ret <- parse$parse_gff_tree(gff_file)
        names(ret) <-
            c(
                "chr_to_gene",
                "transcript_dict",
                "gene_to_transcript",
                "transcript_to_exon"
            )

        ret
    }, args = args)

    ret
}
