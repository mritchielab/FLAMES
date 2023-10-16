#' cutadapt wrapper
#' @description trim TSO adaptor with cutadapt
#' @param args arguments to be passed to cutadapt
#' @return Exit code of cutadapt
#'
#' @examples
#' \dontrun{
#'  cutadapt("-h")
#' }
#' @importFrom reticulate import_from_path
#' @importFrom basilisk basiliskRun
#' @export
cutadapt <- function(args) {
  basiliskRun(env = flames_env, fun = function(x) {
    subprocess <- reticulate::import("subprocess")
    builtin <- reticulate::import_builtins()
    output <- subprocess$check_output(paste("cutadapt", as.list(x), sep=" "), shell=TRUE)
    output_str <- builtin$str(output, encoding="utf-8")
    builtin$print(output_str)

  }, x = args)
}

# cutadapt usage: cutadapt -a 'CCCATGTACTCTGCGTTGATACCACTGCTT' -o cutadapt.fq  --untrimmed-output untrimmed.fq ../main/1k.out.fq
