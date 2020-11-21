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
    callBasilisk(flames_env, function(args) {
        python_path <- system.file("python", package="FlamesR")
    }, args=args)
}

