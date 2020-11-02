#' Title
#'
#' Desc
#'
#' @param v a vector
#'
#' @importFrom basilisk basiliskStart basiliskStop basiliskRun
#' @importFrom reticulate source_python
#'

testFunc <- function(v) {
  proc <- basiliskStart(test_env)
  on.exit(basiliskStop(proc))

  python_path <- system.file("python", package="FlamesR")
  print(python_path)

  tester <- basiliskRun(proc, function(vec) {
    # the R function which can call python functions
    source_python(file=paste(python_path, "test_python.py", sep="/"))
    x <- python_test(vec)

    x
  }, vec = v)

  tester

}
