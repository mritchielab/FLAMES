#' BLAZE Assign reads to cell barcodes. 
#'
#' @description
#' Uses BLAZE to assign reads to cell barcodes. 
#' Uses default options for BLAZE, see BLAZE documentation for details (https://github.com/shimlab/BLAZE).
#'
#' @param config Parsed list of FLAMES config file
#' @param fq_in File path to the fastq file used as a query sequence file
#' @param outdir Output folder
#' @param prefix String, the prefix (e.g. sample name) for the output files
#' @param threads Integer, threads for BLAZE to use, FLAMES will try to detect cores if this parameter is not provided.
#' @param blaze_config List, additional BLAZE configuration parameters
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
#' fastq1_url <- 'https://raw.githubusercontent.com/OliverVoogd/FLAMESData/master/data/fastq/sample1.fastq.gz'
#' fastq1 <- bfc[[names(BiocFileCache::bfcadd(bfc, 'Fastq1', fastq1_url))]]
#' outdir <- tempfile()
#' dir.create(outdir)
#' config = jsonlite::fromJSON(system.file('extdata/blaze_flames.json', package = 'FLAMES'))
#' config$blaze_parameters['output-prefix'] <- outdir
#' blaze(config$blaze_parameters, fastq1)
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
        bc_list_10x <- bfc[[names(BiocFileCache::bfcadd(bfc, 'bc_list_10x', bc_list_10x_url))]]
        blaze_argv <- paste(blaze_argv, '--full-bc-whitelist', bc_list_10x)

        blaze_argv <- paste(blaze_argv, fq_in)

        ret <-
            callBasilisk(blaze_env, function(blaze_argv) {

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
        ret # return the filename of demultiplexed fastq
    }

            #  <fastq directory>
            
            # Required argument:
            #     --expect-cells

            # Options:
            #     -h, --help
            #     --output_fastq
            #     --kit-version <v2 or v3>:
            #     --minQ <INT>:
            #     --threads <INT>
            #     --batch-size <INT>
            #     --full-bc-whitelist <path to file>
            #     --out-putative-bc <filename_prefix>
            #     --out-bc-whitelist <filename_prefix>


            # High sensitivity mode:

            #     --high-sensitivity-mode:
            #         Note that --emptydrop is recommanded specified with this mode (See details below).

            # Empty droplet BCs
            #     --emptydrop
            #     --emptydrop-max-count <INT>