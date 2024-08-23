FROM bioconductor/bioconductor_docker:devel

WORKDIR /home/rstudio

COPY --chown=rstudio:rstudio . /home/rstudio/

RUN Rscript -e "options(repos = c(CRAN = 'https://cran.r-project.org')); BiocManager::install(ask=FALSE)"

RUN Rscript -e "options(repos = c(CRAN = 'https://cran.r-project.org')); devtools::install('.', dependencies=TRUE, build_vignettes=TRUE, repos = BiocManager::repositories())"

RUN sudo apt-get update && sudo apt-get install -y samtools minimap2 && wget -O- https://github.com/attractivechaos/k8/releases/download/v1.2/k8-1.2.tar.bz2 | tar -jx

RUN sudo ln -s /home/rstudio/k8-1.2/k8-x86_64-Linux /bin/k8

USER rstudio

RUN Rscript -e "basilisk::basiliskRun(env = FLAMES:::flames_env, fun = function(){})"

USER root

RUN chmod -R 777 /home/rstudio
