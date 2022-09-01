#' @importFrom basilisk BasiliskEnvironment
flames_env <- BasiliskEnvironment(
    envname = "flames_env", pkgname = "FLAMES",
    packages = c(
        "python>=3.7",
        # "minimap2==2.17",
        "numpy==1.16.5",
        "editdistance==0.5.3",
        "scipy==1.2.0",
        "pysam==0.18.0"
    ),
    channels = c("bioconda", "conda-forge")
)