---
title: "FLAMES"
author: "Oliver Voogd"
package: FLAMES                
output: 
  BiocStyle::html_document
bibliography: ref.bib
vignette: >
  %\VignetteIndexEntry{FLAMES}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}

---


# FLAMES

The FLAMES package provides a framework for performing single-cell and bulk read full-length analysis of mutations and splicing. FLAMES performs cell barcode and UMI assignment from nanopore reads as well as semi-supervised isoform detection and quantification. FLAMES is designed to be an easy and quick to use, powerful workflow for isoform detection and quantification, splicing analysis and mutation detection.

#### image of workflow

Input to FLAMES are fastq files generated from the long-read platform. Using the cell barcode annotation obtained from short-read data as the reference, it identifies and trims cell barcodes/UMI sequences from the long reads. After barcode assignment, all reads were aligned to the relevant genome to obtain a draft read alignment. The draft alignment is then polished and grouped to generate a consensus transcript assembly. All reads are aligned again using the transcript assembly as the reference and quantified. 
### More 'motivating' commentary on the process of the pipeline

This vignette will detail the process of running the FLAMES pipeline. It uses the bulk data pipeline (`bulk_long_pipeline`) as an example, however this vignette will be useful for those running the other FLAMES entry point, the single cell pipeline (`sc_long_pipeline`).

FLAMES uses minimap2 for both the read alignment and read realignment to isoform steps of the pipeline (). If a user wishes to use the FLAMES pipeline without access to minimap2, please refer to the vignette `Vignette for FLAMES on Windows`, which provides instructions for running the pipeline on Windows, which is applicable to any system without minimap2 installed. 

# FLAMES Execution

### Environment setup
To get started, the pipeline needs access to a gene annotation file in gff3 or gtf format, a directory containing one or more fastq files (which will be merged as pre-processing), a genome fasta file, as well as the file path to minimap2, and the file path to the directory to hold output files.

For this vignette, these files are downloaded from GitHub using BiocFileCache.

```r
# download required files using BiocFileCache
temp_path <- tempfile()
bfc <- BiocFileCache::BiocFileCache(temp_path, ask=FALSE)
file_url <- 
  "https://raw.githubusercontent.com/OliverVoogd/FLAMESData/master/data"
annot <- bfc[[names(BiocFileCache::bfcadd(bfc, "Annotation", 
                                          paste(file_url, 
                                                "SIRV_isoforms_multi-fasta-annotation_C_170612a.gtf", 
                                                sep="/")))]] # [[ notation is used to get the local file path of the downloaded file
genome_fa <- bfc[[names(BiocFileCache::bfcadd(bfc, 
                                              "Genomefa", 
                                              paste(file_url, 
                                                    "SIRV_isoforms_multi-fasta_170612a.fasta", 
                                                    sep="/")))]]

# download the two fastq files, move them to a folder to be merged together
fastq1 <- bfc[[names(BiocFileCache::bfcadd(bfc, "Fastq1", paste(file_url, "fastq/sample1.fastq.gz", sep="/")))]]
fastq2 <- bfc[[names(BiocFileCache::bfcadd(bfc, "Fastq2", paste(file_url, "fastq/sample2.fastq.gz", sep="/")))]]
fastq_dir <- paste(temp_path, "fastq_dir", sep="/") # the downloaded fastq files need to be in a directory to be merged together
dir.create(fastq_dir)
file.copy(c(fastq1, fastq2), fastq_dir)
#> [1] TRUE TRUE
unlink(c(fastq1, fastq2)) # the original files can be deleted

# setup other environment variables
config_file <- system.file("extdata/SIRV_config_default.json", package="FLAMES") # the configuration file is included with the FLAMES package
outdir <- tempfile() # create a temporary output directory
if (!dir.exists(outdir)) {
    dir.create(outdir)
}
```

`config_file` is an optional argument which can be given to both `bulk_long_pipeline` and `sc_long_pipeline` in order to customise the execution of the pipelines. It is expected as a JSON file, and an example can be found at /Users/voogd.o/Library/R/4.0/library/FLAMES/extdata/SIRV_config_default.json. 
If `config_file` is not given, the pipeline can instead be customised by altering any of the optional arguments the pipeline allows. These customisations are then stored in a JSON file saved in the specified out directory, which allows for easier reproduction of pipeline execution. More information on the optional arguments at the contents of the JSON configuration file can be found by running `?create_config` and `?bulk_long_pipeline`.
This vignette uses the default configuration file.

### FLAMES execution
Once the environment has been setup, the pipeline can be executed by running:

```r
library(FLAMES)
summarizedExperiment <- bulk_long_pipeline(annot=annot, fastq=fastq_dir, outdir=temp_path, 
                                           genome_fa=genome_fa, minimap2_dir=mm2_dir,
                                           config_file = config_file)
```

### FLAMES termination
The directory `outdir` now contains several output files returned from this pipeline. The output files generated by this pipeline are:

\itemize{
  \item{transcript_count.csv.gz}{ - a transcript count matrix (also contained in the SummarizedExperiment)}
  \item{isoform_annotated.filtered.gff3}{ - isoforms in gff3 format (also contained in the SummarizedExperiment)}
  \item{transcript_assembly.fa}{ - transcript sequence from the isoforms}
  \item{align2genome.bam}{ - sorted BAM file with reads aligned to genome}
  \item{realign2transcript.bam}{ - sorted realigned BAM file using the transcript_assembly.fa as reference}
  \item{tss_tes.bedgraph}{ - TSS TES enrichment for all reads (for QC)}
 }

The pipeline also returns a SummarizedExperiment  or SingleCellExperiment object, depending on the pipeline run, containing the data from `transcript_count.csv.gz`and `isoform_annotated.filtered.gff3`.


## FLAMES on Windows
Due to FLAMES requiring minimap2 to complete the pipeline, the straight FLAMES pipeline functions `bulk_long_pipeline` and `sc_long_pipeline` won't run on a Windows OS. To overcome this issue, the vignette 'Vignette for FLAMES on Windows' describes the alternate method of running the FLAMES pipelines, which requires acccess to minimap2 on an external system

# Session Info

```
#> R version 4.0.5 (2021-03-31)
#> Platform: x86_64-apple-darwin17.0 (64-bit)
#> Running under: macOS Catalina 10.15.7
#> 
#> Matrix products: default
#> BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
#> LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib
#> 
#> locale:
#> [1] en_AU.UTF-8/en_AU.UTF-8/en_AU.UTF-8/C/en_AU.UTF-8/en_AU.UTF-8
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> loaded via a namespace (and not attached):
#>  [1] Rcpp_1.0.6           compiler_4.0.5       pillar_1.6.1        
#>  [4] dbplyr_2.1.1         tools_4.0.5          digest_0.6.27       
#>  [7] bit_4.0.4            BiocFileCache_1.14.0 evaluate_0.14       
#> [10] RSQLite_2.2.7        memoise_2.0.0        lifecycle_1.0.0     
#> [13] tibble_3.1.2         pkgconfig_2.0.3      rlang_0.4.11        
#> [16] DBI_1.1.1            curl_4.3.1           xfun_0.22           
#> [19] fastmap_1.1.0        dplyr_1.0.6          stringr_1.4.0       
#> [22] httr_1.4.2           knitr_1.32           generics_0.1.0      
#> [25] vctrs_0.3.8          rappdirs_0.3.3       bit64_4.0.5         
#> [28] tidyselect_1.1.1     glue_1.4.2           R6_2.5.0            
#> [31] fansi_0.4.2          rmarkdown_2.7        purrr_0.3.4         
#> [34] blob_1.2.1           magrittr_2.0.1       htmltools_0.5.1.1   
#> [37] ellipsis_0.3.2       assertthat_0.2.1     utf8_1.2.1          
#> [40] stringi_1.6.1        cachem_1.0.5         crayon_1.4.1
```

# References
