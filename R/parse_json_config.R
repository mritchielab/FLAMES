#' Title
#'
#' DESC
#' 
#' @param name desc
#'
#' @param name desc
#' @importFrom reticulate import_from_path
#' @export
parse_json_config <- function(json_config) {
    config <- callBasilisk(flames_env, function(config) {
        python_path <- system.file("python", package="FlamesR")
        js <- reticulate::import_from_path("parse_config", python_path)
        js$parse_json_config(config)
    }, config=json_config)

    config
}

