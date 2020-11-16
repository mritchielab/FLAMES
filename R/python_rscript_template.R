#' Title
#'
#' DESC
#' 
#' @param name desc
#'
#' @param name desc
#' @importFrom reticulate import_from_path
#' @export
name <- function(args) {
    python_path <- system.file("python", package="FlamesR")
    callBasilisk(flames_env, FUN)
}

