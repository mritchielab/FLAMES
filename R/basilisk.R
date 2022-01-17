#' @importFrom basilisk BasiliskEnvironment
flames_env <- BasiliskEnvironment(
    envname = "flames_env", pkgname = "FLAMES",
    packages = c(
        "python==2.7.15.0",
        # "minimap2==2.17",
        "numpy==1.16.5",
        "editdistance==0.5.3",
        "bamnostic==1.1.7"
    ),
    channels = c("bioconda", "conda-forge")
)