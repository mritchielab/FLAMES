# flames_nopysam_env <- BasiliskEnvironment(
#    envname = "flames_nopysam_env", pkgname = "FLAMES",
#    packages = c(
#        "python==2.7.15.0",
#        # "minimap2==2.17",
#        "numpy==1.16.5",
#        "editdistance==0.5.3",
#        "bamnostic==1.1.7"
#    ),
#    channels = c("bioconda", "conda-forge")
# )

#' @importFrom basilisk BasiliskEnvironment
flames_env <- BasiliskEnvironment(
    envname = "flames_env", pkgname = "FLAMES",
    packages = c(
        "python==3.7",
        # "minimap2==2.17",
        "numpy==1.16.5",
        "editdistance==0.5.3",
        "scipy==1.2.0",
        "pysam==0.18.0"
))

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
# flames_env <- BasiliskEnvironment(envname="full_env",
#     pkgname="FLAMES",
#     packages=c("python==2.7.15.0",
#         "pysam==0.16.0.1",
#         "minimap2==2.17",
#         "numpy==1.16.5",
#         "editdistance==0.5.3",
#         "bzip2==1.0.8",
#         "c-ares==1.11.0",
#         "ca-certificates==2020.11.8",
#         "certifi==2019.11.28",
#         "htslib==1.11",
#         "k8==0.2.5",
#         "krb5==1.17.1",
#         "libblas==3.9.0",
#         "libcblas==3.9.0",
#         "libcurl==7.71.1",
#         "libcxx==11.0.0",
#         "libdeflate==1.6",
#         "libedit==3.1.20191231",
#         "libev==4.33",
#         "libffi==3.2.1",
#         "libgfortran==3.0.0",
#         "libgfortran5==9.3.0",
#         "liblapack==3.9.0",
#         "libnghttp2==1.41.0",
#         "libopenblas==0.3.12",
#         "libssh2==1.9.0",
#         "llvm-openmp==11.0.0",
#         "ncurses==6.2",
#         "openssl==1.1",
#         #"pip==20.1.1",
#         "python_abi==2.7",
#         "readline==8.0",
#         "setuptools==44.0.0",
#         "sqlite==3.33.0",
#         "tk==8.6.10",
#         "wheel==0.35.1",
#         "xz==5.2.5",
#         "zlib==1.2.11"),
#         channels=c("bioconda", "conda-forge"))
