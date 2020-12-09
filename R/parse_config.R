#' Title
#'
#' DESC
#' 
#' @param name desc
#'
#' @param name desc
#' @importFrom basilisk basiliskStart basiliskStop basiliskRun
#' @importFrom reticulate import_from_path
#' @export
parse_json_config <- function(json_file) {
    callBasilisk(flames_env, function(json) {
        
    }, json=json_file)
}

#' Title
#'
#' DESC
#' 
#' @param name desc
#'
#' @param name desc
#' @importFrom reticulate import_from_path
#' @export
print_config <- function(decoded_dict) {
    callBasilisk(flames_env, function(decoded_dict) {
        python_path <- system.file("python", package="FlamesR")

        conf <- reticulate::import_from_path("parse_config", python_path)
        conf$print_config(decoded_dict)
    }, decoded_dict=decoded_dict)
    invisible()
}

