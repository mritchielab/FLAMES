#' Parse Json Configuration File
#'
#' @description Convert a json configuration file into a named R list, grouped into sub lists according to their 
#'      usage in the Flames pipeline.
#' 
#' @param json_file The file name to convert into an R list.
#'
#' @return A named R list of the parameters in \code{json_file}. Subsections are: \code{pipeline_parameters},
#'      \code{global_parameters}, \code{isoform_parameters}, \code{alignment_parameters}, \code{realign_parameters} and
#'      \code{transcript_counting}.
#' @importFrom basilisk basiliskStart basiliskStop basiliskRun
#' @importFrom reticulate import_from_path
#' @export
parse_json_config <- function(json_file) {
    config <- callBasilisk(flames_env, function(json) {
        python_path <- system.file("python", package="FlamesR")

        conf <- reticulate::import_from_path("parse_config", python_path)

        conf$parse_json_config(json)
    }, json=json_file)

    config
}

#' Print Configuration File
#'
#' @details Print the configuration file, represented as a named list used for the Flames pipeline.
#' 
#' @param config List; the configuration list to print.
#'
#' @importFrom reticulate import_from_path
#' @export
print_config <- function(config) {
    callBasilisk(flames_env, function(config) {
        python_path <- system.file("python", package="FlamesR")

        conf <- reticulate::import_from_path("parse_config", python_path)
        conf$print_config(config)
    }, config=config)
    invisible()
}

#' Write Configuration Dictionary to File
#'
#' @details Print the configuration file, represented as a named list used for the Flames pipeline.
#' 
#' @param config List; the configuration list to print.
#' @param config_file the file to output \code{config} to. Should be .json extension
#' @importFrom reticulate import_from_path
write_config <- function(config, config_file) {
    # write the config file to given file path
    callBasilisk(flames_env, function(config, config_file) {
        python_path <- system.file("python", package="FlamesR")

        conf <- reticulate::import_from_path("parse_config", python_path)
        conf$write_config(config, config_file)
    }, config=config, config_file=config_file)
    invisible()
}

