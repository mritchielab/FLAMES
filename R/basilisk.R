# change the name of this env to something more useful
#' @importFrom basilisk BasiliskEnvironment
flames_env <- BasiliskEnvironment(envname="flames_env",
                                  pkgname="FlamesR",
                                  packages=c("python==2.7", "samtools==1.11", "pysam==0.16.0.1", "minimap2==2.17",
                                            "numpy==1.19.2", "editdistance==0.5.3"),
                                  channels=c("bioconda", "conda-forge"))

test_env <- BasiliskEnvironment(envname="test_env",
                                pkgname="FlamesR",
                                packages=c("python==2.7"))
