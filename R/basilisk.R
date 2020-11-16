# change the name of this env to something more useful
#' @importFrom basilisk BasiliskEnvironment
full_env <- BasiliskEnvironment(envname="full_env",
    pkgname="FlamesR",
    packages=c("python==2.7", "samtools==1.11", "pysam==0.16.0.1", "minimap2==2.17", "numpy==1.16.5", "editdistance==0.5.3",
        "bzip2==1.0.8", "c-ares==1.11.0", "ca-certificates==2020.11.8", "certifi==2019.11.28", "htslib==1.11", "k8==0.2.5",
        "krb5==1.17.1", "libblas==3.9.0", "libclas==3.9.0", "libcurl==7.71.1", "libcxx==11.0.0", "libdeflate==1.6",
        "libedit==3.1.20191231", "libev==4.33", "libffi==3.2.1", "libgfortran==5.0.0", "libgfortran5==9.3.0",
        ))

#' @importFrom basilisk BasiliskEnvironment
#flames_env <- BasiliskEnvironment(envname="flames_env",
#                                  pkgname="FlamesR",
#                                  packages=c("python==2.7", "samtools==1.11", "pysam==0.16.0.1", "minimap2==2.17",
#                                            "numpy==1.16.5", "editdistance==0.5.3", "utils==2.4.0"),
#                                  channels=c("bioconda", "conda-forge"))


#  libcblas           conda-forge/osx-64::libcblas-3.9.0-2_openblas
#  libcurl            conda-forge/osx-64::libcurl-7.71.1-h9bf37e3_8
#  libcxx             conda-forge/osx-64::libcxx-11.0.0-h439d374_0
#  libdeflate         conda-forge/osx-64::libdeflate-1.6-h0b31af3_0
#  libedit            conda-forge/osx-64::libedit-3.1.20191231-h0678c8f_2
#  libev              conda-forge/osx-64::libev-4.33-haf1e3a3_1
#  libffi             bioconda/osx-64::libffi-3.2.1-1
#  libgfortran        conda-forge/osx-64::libgfortran-5.0.0-h7cc5361_13
#  libgfortran5       conda-forge/osx-64::libgfortran5-9.3.0-h7cc5361_13
#  liblapack          conda-forge/osx-64::liblapack-3.9.0-2_openblas
#  libnghttp2         conda-forge/osx-64::libnghttp2-1.41.0-h8a08a2b_1
#  libopenblas        conda-forge/osx-64::libopenblas-0.3.12-openmp_h54245bb_1
#  libssh2            conda-forge/osx-64::libssh2-1.9.0-h8a08a2b_5
#  llvm-openmp        conda-forge/osx-64::llvm-openmp-11.0.0-h73239a0_1
#  minimap2           bioconda/osx-64::minimap2-2.17-hbbe82c9_3
#  ncurses            conda-forge/osx-64::ncurses-6.2-h2e338ed_4
#  numpy              conda-forge/osx-64::numpy-1.16.5-py27hde6bac1_0
#  openssl            conda-forge/osx-64::openssl-1.1.1h-haf1e3a3_0
#  pip                conda-forge/noarch::pip-20.1.1-pyh9f0ad1d_0
#  pysam              bioconda/osx-64::pysam-0.16.0.1-py27hdf4340f_1
#  python             conda-forge/osx-64::python-2.7.15-h8e446fc_1011_cpython
#  python_abi         conda-forge/osx-64::python_abi-2.7-1_cp27m
#  readline           conda-forge/osx-64::readline-8.0-h0678c8f_2
#  samtools           bioconda/osx-64::samtools-1.11-h725deca_0
#  setuptools         conda-forge/osx-64::setuptools-44.0.0-py27_0
#  sqlite             conda-forge/osx-64::sqlite-3.33.0-h960bd1c_1
#  tk                 conda-forge/osx-64::tk-8.6.10-hb0a8c7a_1
#  wheel              conda-forge/noarch::wheel-0.35.1-pyh9f0ad1d_0
#  xz                 conda-forge/osx-64::xz-5.2.5-haf1e3a3_1
#  zlib               conda-forge/osx-64::zlib-1.2.11-h7795811_1010


#' @importFrom basilisk BasiliskEnvironment
#pysam_env <- BasiliskEnvironment(envname="pysam_env", 
#                                 pkgname="FlamesR",
#                                 packages=c("python==2.7", "pysam==0.14", "utils==2.4.0"),
#                                 channels=c("bioconda"))


new_env2 <- BasiliskEnvironment(envname="new_env2",
                               pkgname="FlamesR",
                               packages=c("python==2.7", "pysam==0.15.4", "samtools==1.7"),
                               channels=c("bioconda", "conda-forge"))