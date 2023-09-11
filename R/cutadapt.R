#' cutadapt wrapper
#' @description trim TSO adaptor with cutadapt
#' @param args arguments to be passed to cutadapt
#' @return Exit code of cutadapt
#'
#' @examples
#' \dontrun{
#'  cutadapt("-h")
#' }
#' @export
cutadapt <- function(args) {
  callBasilisk(flames_env, function(x) {
    # python_path <- system.file("python", package = "FLAMES")
    # mod <- reticulate::import_from_path("cutadapt_wrapper", python_path)
    # return(mod$wrapper(as.list(x)))
    subprocess <- reticulate::import("subprocess")
    builtin <- reticulate::import_builtins()
    output <- subprocess$check_output(paste("cutadapt", as.list(x), sep=" "), shell=TRUE)
    output_str <- builtin$str(output, encoding="utf-8")
    builtin$print(output_str)

  }, x = args)
}

#cutadapt -a 'CCCATGTACTCTGCGTTGATACCACTGCTT' -o cutadapt.fq  --untrimmed-output untrimmed.fq ../main/1k.out.fq
