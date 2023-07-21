#' @importFrom basilisk BasiliskEnvironment
flames_env <- BasiliskEnvironment(
    envname = "flames_env", pkgname = "FLAMES",
    packages = c(
        "python==3.10.12",
        # "minimap2==2.17",
        "numpy==1.25.0",
        "editdistance==0.6.2",
        "scipy==1.11.1",
        "pysam==0.21.0",
        "cutadapt==4.4"
    ),
    channels = c("conda-forge", "bioconda")
)
