#' Title
#'
#' DESC
#' 
#' @param name desc
#'
#' @param name desc
#' @export
test <- function() {
    callBasilisk(full_env, function() {
        x <- 10
        python_path <- system.file("python", package="FlamesR")

        a <- reticulate::import_from_path("test2", python_path)
        #print(a$add(10, 20))
        # Thank fuck importing pysam finally works!!
        a$add(10, 20)
    })
}

