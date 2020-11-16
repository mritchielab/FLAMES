#' Title
#'
#' DESC
#' 
#' @param minimap2_prog_path Directory containing minimap2, k8 and paftools.js. These are 
#' used to convert gff3 to bed12.
#'
#' @param gff3_file The gff3_file to convert
#' @param bed12_file The filename of the bed12 output file.
#' @importFrom reticulate import_from_path source_python
#' @export
gff3_to_bed12 <- function(minimap2_prog_path, gff3_file, bed12_file) {
    python_path <- system.file("python", package="FlamesR")
    print(python_path)
    callBasilisk(flames_env, function(mm2_path, gff3, bed12) {
        reticulate::import_from_path("test2", python_path)
        reticulate::source_python(paste(python_path, "test2.py", sep=.Platform$file.sep))
        
        add(10, 20)
    #    align <-reticulate::import_from_path("minimap2_align", python_path)
     #   reticulate::source_python(paste(python_path, "minimap2_align.py", sep=.Platform$file.sep))

        #align$gff3_to_bed12(mm2_path, gff3, bed12)
    }, mm2_path=minimap2_prog_path, gff3=gff3_file, bed12=bed12_file)

    10
}


#' Title
#'
#' DESC
#' 
#' @param name desc
#'
#' @param name desc
#' @export
test <- function() {
    callBasilisk(new_env2, function() {
        
        10
    })
}

