#' @importFrom basilisk BasiliskEnvironment
flames_env <- BasiliskEnvironment(
    envname = "flames_env", pkgname = "FLAMES",
    pip = c("fast-edit-distance==1.2.1", "blaze2==2.2.*", "matplotlib==3.5.3"),
    packages = c(
        "python==3.10",
        "numpy==1.25.0",
        "scipy==1.11.1",
        "pysam==0.21.0",
        "cutadapt==4.4",
        "tqdm==4.64.1",
        "pandas==1.3.5"
    ),
    channels = c("conda-forge", "bioconda", "defaults")
)