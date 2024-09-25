#' @importFrom basilisk BasiliskEnvironment
flames_env <- BasiliskEnvironment(
    envname = "flames_env", pkgname = "FLAMES",
    pip = c("fast-edit-distance==1.2.1", "blaze2==2.2.*", "matplotlib==3.9.1"),
    packages = c(
        "python==3.11.9",
        "numpy==2.1.1",
        "scipy==1.14.1",
        "pysam==0.22.1",
        "cutadapt==4.9",
        "tqdm==4.66.5",
        "pandas==2.2.3"
    ),
    channels = c("conda-forge", "bioconda")
)

#' @importFrom basilisk BasiliskEnvironment
bins_env <- BasiliskEnvironment(
    envname = "bins_env", pkgname = "FLAMES",
    packages = c(
        "python==3.8.20",
        "oarfish==0.6.2",
        "minimap2==2.28",
        "samtools==1.21",
        "k8==1.2"
    ),
    channels = c("conda-forge", "bioconda")
)
