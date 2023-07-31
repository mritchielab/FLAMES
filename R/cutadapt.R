#' cutadapt wrapper
#' @description trim TSO adaptor with cutadapt
#' @param args arguments to be passed to cutadapt
#' @return Exit code of cutadapt
#'
#' @examples
#' cutadapt("-h")
#'
#' @export
cutadapt <- function(args) {
  callBasilisk(FLAMES:::flames_env, function(x) {
    python_path <- system.file("python", package = "FLAMES")
    mod <- reticulate::import_from_path("cutadapt_wrapper", python_path)
    return(mod$wrapper(as.list(x)))
  }, x = args)
}

#cutadapt -a 'CCCATGTACTCTGCGTTGATACCACTGCTT' -o cutadapt.fq  --untrimmed-output untrimmed.fq ../main/1k.out.fq
