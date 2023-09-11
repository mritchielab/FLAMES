#' BLAZE Assign reads to cell barcodes. 
#'
#' @description
#' Uses BLAZE to generate barcode list and assign reads to cell barcodes. 
#' Uses default options for BLAZE, see BLAZE documentation for details (https://github.com/shimlab/BLAZE).
#'
#' @param blaze_config List, additional BLAZE configuration parameters
#' @param fq_in File path to the fastq file used as a query sequence file
#'
#' @return a \code{data.frame} summarising the reads aligned
#'
#' @importFrom parallel detectCores
#' @export
#' @examples
#' temp_path <- tempfile()
#' bfc <- BiocFileCache::BiocFileCache(temp_path, ask = FALSE)
#' bc_list_10x_url <- 'https://github.com/shimlab/BLAZE/blob/main/10X_bc/3M-february-2018.zip'
#' bc_list_10x <- bfc[[names(BiocFileCache::bfcadd(bfc, 'bc_list_10x', bc_list_10x_url))]]
#' fastq1_url <- 'https://raw.githubusercontent.com/shimlab/BLAZE/main/test/data/FAR20033_pass_51e510db_100.fastq'
#' fastq1 <- bfc[[names(BiocFileCache::bfcadd(bfc, 'Fastq1', fastq1_url))]]
#' outdir <- tempfile()
#' dir.create(outdir)
#' config = jsonlite::fromJSON(system.file('extdata/blaze_flames.json', package = 'FLAMES'))
#' config$blaze_parameters['output-prefix'] <- outdir
#' \dontrun{
#'    blaze(config$blaze_parameters, fastq1)
#' }
#' @importFrom reticulate import_from_path dict
#' @export
blaze <- function(blaze_config, fq_in) {
        
        # command line arguments for blaze
        blaze_argv <- paste("")

        if (blaze_config['overwrite'] == TRUE) {
            blaze_argv <- paste(blaze_argv, '--overwrite ')
        }
        blaze_config['overwrite'] <- NULL
        for (arg in names(blaze_config)) {
            blaze_argv <- paste(blaze_argv, paste0('--',arg), blaze_config[arg])}
        
        # prepare 10X whitelist
        temp_path <- tempfile()
        bfc <- BiocFileCache::BiocFileCache(temp_path, ask = FALSE)
        bc_list_10x_url <- 'https://github.com/shimlab/BLAZE/raw/main/10X_bc/3M-february-2018.zip'
        cat('Downloading the full whitelist from 10X...')
        bc_list_10x <- bfc[[names(BiocFileCache::bfcadd(x=bfc, rname='bc_list_10x', fpath=bc_list_10x_url))]]
        blaze_argv <- paste(blaze_argv, '--full-bc-whitelist', bc_list_10x)

        blaze_argv <- paste(blaze_argv, fq_in)

        ret <-
            callBasilisk(flames_env, function(blaze_argv) {

                blaze_path <- system.file("blaze", package = "FLAMES")
                cat("Running BLAZE...\n")
                cat("Argument: ", blaze_argv, "\n")
                blaze <-
                    reticulate::import_from_path("blaze", blaze_path)
                ret <-
                    blaze$main(blaze_argv)

                ret
            }, blaze_argv = blaze_argv
            )
        #ret # return the filename of demultiplexed fastq
    }
