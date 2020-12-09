#' Title
#'
#' DESC
#' 
#' @param name desc
#'
#' @param name desc
#' @import reticulate import_from_path py_dict
#' @export
test <- function(first, ...) {
    kwargs=reticulate::dict(...)
    callBasilisk(flames_env, function(first, kwargs) {
        python_path <- system.file("python", package="FlamesR")
        print(kwargs)
        a <- reticulate::import_from_path("test2", python_path)
        a$t3(first, kwargs, a=12, b=54, c="hello")
        #print(chr_to_gene)
    }, first=first, kwargs=kwargs)
}
#' Title
#'
#' DESC
#' 
#' @param name desc
#'
#' @param name desc
#' @export
test2 <- function(g) {
    callBasilisk(flames_env, function(g) {
        python_path <- system.file("python", package="FlamesR")

        a <- reticulate::import_from_path("test2", python_path)
        juncs <- a$t(g)
        a$t2(juncs)
        10
    }, g=g)
}

test_parse <- function(annot) {
    callBasilisk(flames_env, function(a) {
        p <- system.file("python", package="FlamesR")
        c <- reticulate::import_from_path("parse_gene_anno", p)
        c$parse_gff_tree(a)
    }, a=annot)
}

