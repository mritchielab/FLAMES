#' Default Configuration File
#'
#' @description file path to the default FLAMES configuration file
#' @examples
#' config <- get_default_config_file()
#' @return file path to the FLAMES default configuration file
#' @export
get_default_config_file <- function() {
    system.file("extdata/SIRV_config_default.json", package = "FLAMES")
}
