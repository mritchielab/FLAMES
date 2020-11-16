#' Title
#'
#' Desc
#'
#' @param v a vector
#'
#' @importFrom basilisk basiliskStart basiliskStop basiliskRun
#' @importFrom reticulate source_python
testFunc <- function() {
  proc <- basiliskStart(pysam_env)
  on.exit(basiliskStop(proc))

  python_path <- system.file("python", package="FlamesR")
  print(python_path)

  tester <- basiliskRun(proc, function(x, y) {
    # the R function which can call python functions
    #reticulate::source_python((file=paste(python_path, "test2.py", sep="/")))
    reticulate::import_from_path("test_python", python_path)
    reticulate::source_python(file=paste(python_path, "test_python.py", sep="/"))
    x <- python_test2(x, y)

    x
  }, x = 1, y = 2)

  #tester

}

#' Title
#'
#' Desc
#'
#' @param v a vector
#'
#' @importFrom basilisk basiliskStart basiliskStop basiliskRun
#' @importFrom reticulate source_python
#'
tRun <- function(proc) {
  t <- basiliskRun(proc, function() {
    source_python(paste(system.file("inst/python", package="FlamesR"), "test_python.py", sep="/"))
    python_test(c(1,2,3,4,5))
  })

  t
}
