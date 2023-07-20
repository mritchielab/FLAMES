#' @importFrom basilisk BasiliskEnvironment
flames_env <- BasiliskEnvironment(
    envname = "flames_env", pkgname = "FLAMES",
    packages = c(
        "python==3.10.12",
        # "minimap2==2.17",
        "numpy==1.25.0",
        "editdistance==0.6.2",
        "scipy==1.11.1",
        "pysam==0.21.0"
    ),
    channels = c("conda-forge", "bioconda")
)

blaze_env <- BasiliskEnvironment(
    envname = "blaze_env", pkgname = "FLAMES",
    pip = c("fast-edit-distance==1.2.0"),
    channels = c('conda-forge','bioconda', 'defaults'),
    packages = c(
        "python==3.7",
        "biopython==1.79", #blaze specific
        "pandas==1.3.5",#blaze specific
        "numpy==1.21.6", #blaze specific
        "matplotlib==3.5.3", #blaze specific
        "tqdm==4.64.1"

))
