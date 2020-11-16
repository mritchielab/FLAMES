#' Title
#'
#' Semi-supervised isofrom detection and annotation from long read data.
#' This variant is meant for bulk samples. THIS FUNCTION IS GOING TO BE THE ENTIRE 
#' FLAMES PROCESS. IT SHOULD BE BUILT FROM THE SMALLER FUNCTIONS THAT PREVIOUSLY WERE
#' CARRIED OUT UNDER THE HOOD. MOST FUNCTIONS THAT THIS ONE USES SHOULD BE AVAILABLE TO THE USER.
#'
#'  THIS NEED PROPER TYPES
#' @param annot gene annotations in gff3 file format. str
#'
#' @param fastq_dir the folder containing fastq files, each containing data from one sample. str
#'
#' @param bam aligned bam file (sorted and indexed). Overwrite --infq? str
#'
#' @param outdir directory to deposite all results in rootdir. Use absolute path? str
#'
#' @param genome_fa genome fasta file. str
#'
#' @param minimap2_dir directory containing minimap2, k8 and paftools.js program. k8 and paftools.js are used to convert gff3 to bed12. str
#'
#' @param config_file json configuration files. str
#'
#' @param downsample_ratio downsampling ratio if performing downsampling analysis. str
#' @importFrom reticulate import_from_path
#' @export
bulk_long_pipeline <- function(annot, fastq_dir, bam, outdir, genome_fa,
                                minimap2_dir, config_file, downsample_ratio) {
    # most of this is transcribed from bulk_long_pipeline
    # argument validation
    infq <- paste(outdir, "merged.fastq.gz", sep="/")
    bc_file <- paste(outdir, "pseudo_barcode_annotation.csv", sep="/")

    print("Preprocessing bulk fastqs...")
    # create output directory if one doesn't exist
    if (!dir.exists(outdir)) {
        print("Output directory does not exists: one is being created")
        dir.create(outdir)
    }

    # run the merge_bulk_fastq function as preprocessing
    merge_bulk_fastq(fastq_dir, bc_file, infq)

    # basilisk setup. flames_env contains all modules needed to run FLAMES python code
    callBasilisk(flames_env, function(a, i, b, o, fa, m, c,
                                        d, bc_f, infastq) {
        print('here')
        python_path <- system.file("python", package="FlamesR")

        bulk_long <-
            reticulate::import_from_path("bulk_long_pipeline", python_path)

        #source_python(file=paste(py, "bulk_long_pipeline.py", sep="/"))
        print("Running FLAMES pipeline...")
        # call the python function. DO THESE ARGUMENTS NEED CONVERSION?
        bulk_long$bulk_long_pipeline(a, i, b, o, f, m, c, d, bc_f, infastq)

        print("finished")
        10
    }, a=annot, i=fastq_dir, b=bam, o=outdir, fa=genome_fa, m=minimap2_dir,
        c=config_file, d=downsample_ratio, bc_f=bc_file, infastq=infq)
}