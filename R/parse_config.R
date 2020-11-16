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

