# Package index

## About FLAMES

- [`FLAMES`](https://mritchielab.github.io/FLAMES/reference/FLAMES.md) :
  FLAMES: full-length analysis of mutations and splicing

## Pipelines

- [`BulkPipeline()`](https://mritchielab.github.io/FLAMES/reference/BulkPipeline.md)
  : Pipeline for bulk long read RNA-seq data processing
- [`SingleCellPipeline()`](https://mritchielab.github.io/FLAMES/reference/SingleCellPipeline.md)
  : Pipeline for Single Cell Data
- [`MultiSampleSCPipeline()`](https://mritchielab.github.io/FLAMES/reference/MultiSampleSCPipeline.md)
  : Pipeline for multi-sample long-read scRNA-seq data
- [`example_pipeline()`](https://mritchielab.github.io/FLAMES/reference/example_pipeline.md)
  : Example pipelins

## Pipeline execution

- [`run_FLAMES()`](https://mritchielab.github.io/FLAMES/reference/run_FLAMES.md)
  : Execute a FLAMES pipeline
- [`resume_FLAMES()`](https://mritchielab.github.io/FLAMES/reference/resume_FLAMES.md)
  : Resume a FLAMES pipeline
- [`run_step()`](https://mritchielab.github.io/FLAMES/reference/run_step.md)
  : Execute a single step of the FLAMES pipeline

## Inspecting pipelines

- [`steps()`](https://mritchielab.github.io/FLAMES/reference/steps.md) :
  Steps to perform in the pipeline
- [`` `steps<-`() ``](https://mritchielab.github.io/FLAMES/reference/steps-set.md)
  : Set steps to perform in the pipeline
- [`config()`](https://mritchielab.github.io/FLAMES/reference/config.md)
  : Get pipeline configurations
- [`` `config<-`() ``](https://mritchielab.github.io/FLAMES/reference/config-set.md)
  : Set pipeline configurations
- [`controllers()`](https://mritchielab.github.io/FLAMES/reference/controllers.md)
  : Get controllers
- [`` `controllers<-`() ``](https://mritchielab.github.io/FLAMES/reference/controllers-set.md)
  : Set controllers
- [`experiment()`](https://mritchielab.github.io/FLAMES/reference/experiment.md)
  : Get pipeline results
- [`plot_durations()`](https://mritchielab.github.io/FLAMES/reference/plot_durations.md)
  : Plot pipeline step durations

## Barcode demultiplexing

- [`find_barcode()`](https://mritchielab.github.io/FLAMES/reference/find_barcode.md)
  : Match Cell Barcodes
- [`barcode_segment()`](https://mritchielab.github.io/FLAMES/reference/barcode_segment.md)
  : Create a Flexiplex barcode segment
- [`barcode_group()`](https://mritchielab.github.io/FLAMES/reference/barcode_group.md)
  : Create a Flexiplex barcode group
- [`blaze()`](https://mritchielab.github.io/FLAMES/reference/blaze.md) :
  BLAZE Assign reads to cell barcodes.
- [`flexiplex()`](https://mritchielab.github.io/FLAMES/reference/flexiplex.md)
  : Rcpp port of flexiplex
- [`cutadapt()`](https://mritchielab.github.io/FLAMES/reference/cutadapt.md)
  : cutadapt wrapper

## Alignment

- [`annotation_to_fasta()`](https://mritchielab.github.io/FLAMES/reference/annotation_to_fasta.md)
  : GTF/GFF to FASTA conversion
- [`index_genome()`](https://mritchielab.github.io/FLAMES/reference/index_genome.md)
  : Index the reference genome for minimap2

## Transcript identification

- [`find_isoform()`](https://mritchielab.github.io/FLAMES/reference/find_isoform.md)
  : Isoform identification

## Quantification

- [`quantify_gene()`](https://mritchielab.github.io/FLAMES/reference/quantify_gene.md)
  : Gene quantification
- [`quantify_transcript()`](https://mritchielab.github.io/FLAMES/reference/quantify_transcript.md)
  : Transcript quantification
- [`quantify_transcript_flames()`](https://mritchielab.github.io/FLAMES/reference/quantify_transcript_flames.md)
  : FLAMES Transcript quantification

## Mutation analysis

- [`find_variants()`](https://mritchielab.github.io/FLAMES/reference/find_variants.md)
  : bulk variant identification
- [`sc_mutations()`](https://mritchielab.github.io/FLAMES/reference/sc_mutations.md)
  : Variant count for single-cell data
- [`sc_genotype()`](https://mritchielab.github.io/FLAMES/reference/sc_genotype.md)
  : Genotype a single-cell mutation
- [`mutation_positions()`](https://mritchielab.github.io/FLAMES/reference/mutation_positions.md)
  : Calculate mutation positions within the gene body

## Visualization

- [`plot_coverage()`](https://mritchielab.github.io/FLAMES/reference/plot_coverage.md)
  : plot read coverages
- [`plot_demultiplex()`](https://mritchielab.github.io/FLAMES/reference/plot_demultiplex.md)
  : Plot Cell Barcode demultiplex statistics
- [`plot_isoform_heatmap()`](https://mritchielab.github.io/FLAMES/reference/plot_isoform_heatmap.md)
  : FLAMES heetmap plots
- [`plot_isoform_reduced_dim()`](https://mritchielab.github.io/FLAMES/reference/plot_isoform_reduced_dim.md)
  : FLAMES isoform reduced dimensions plots
- [`plot_isoforms()`](https://mritchielab.github.io/FLAMES/reference/plot_isoforms.md)
  : Plot isoforms
- [`sc_plot_genotype()`](https://mritchielab.github.io/FLAMES/reference/sc_plot_genotype.md)
  : Plot genotype of single-cell data

## Alignment coverages

- [`get_coverage()`](https://mritchielab.github.io/FLAMES/reference/get_coverage.md)
  : Get read coverages from BAM file
- [`filter_coverage()`](https://mritchielab.github.io/FLAMES/reference/filter_coverage.md)
  : Filter transcript coverage
- [`convolution_filter()`](https://mritchielab.github.io/FLAMES/reference/convolution_filter.md)
  : Convolution filter for smoothing transcript coverages
- [`weight_transcripts()`](https://mritchielab.github.io/FLAMES/reference/weight_transcripts.md)
  : Weight transcripts by read counts

## Analysis of FLT-seq sub-samlped data

- [`combine_sce()`](https://mritchielab.github.io/FLAMES/reference/combine_sce.md)
  : Combine SCE
- [`sc_impute_transcript()`](https://mritchielab.github.io/FLAMES/reference/sc_impute_transcript.md)
  : Impute missing transcript counts

## Analysis of spatial transcriptomics data

- [`create_spe()`](https://mritchielab.github.io/FLAMES/reference/create_spe.md)
  : Create a SpatialExperiment object
- [`plot_spatial_isoform()`](https://mritchielab.github.io/FLAMES/reference/plot_spatial_isoform.md)
  : Plot spatial pie chart of isoforms
- [`plot_spatial_pie()`](https://mritchielab.github.io/FLAMES/reference/plot_spatial_pie.md)
  : Plot spatial pie chart
- [`plot_spatial_feature()`](https://mritchielab.github.io/FLAMES/reference/plot_spatial_feature.md)
  : Plot feature on spatial image

## Analysis of single-cell data

- [`find_diversity()`](https://mritchielab.github.io/FLAMES/reference/find_diversity.md)
  : Compute Gene Isoform Entropy Matrix
- [`sc_DTU_analysis()`](https://mritchielab.github.io/FLAMES/reference/sc_DTU_analysis.md)
  : FLAMES Differential Transcript Usage Analysis

## Miscellaneous

- [`add_gene_counts()`](https://mritchielab.github.io/FLAMES/reference/add_gene_counts.md)
  :

  Add gene counts to a `SingleCellExperiment` object

- [`create_config()`](https://mritchielab.github.io/FLAMES/reference/create_config.md)
  : Create Configuration File From Arguments

- [`load_config()`](https://mritchielab.github.io/FLAMES/reference/load_config.md)
  : Load Configurations

- [`create_sce_from_dir()`](https://mritchielab.github.io/FLAMES/reference/create_sce_from_dir.md)
  :

  Create `SingleCellExperiment` object from `FLAMES` output folder

- [`create_se_from_dir()`](https://mritchielab.github.io/FLAMES/reference/create_se_from_dir.md)
  :

  Create `SummarizedExperiment` object from `FLAMES` output folder

- [`filter_annotation()`](https://mritchielab.github.io/FLAMES/reference/filter_annotation.md)
  : filter annotation for plotting coverages

- [`get_GRangesList()`](https://mritchielab.github.io/FLAMES/reference/get_GRangesList.md)
  : Parse FLAMES' GFF output

- [`demultiplex_sockeye()`](https://mritchielab.github.io/FLAMES/reference/demultiplex_sockeye.md)
  : Demultiplex reads using Sockeye outputs

- [`scmixology_lib10`](https://mritchielab.github.io/FLAMES/reference/scmixology_lib10.md)
  : scMixology short-read gene counts - sample 2

- [`scmixology_lib10_transcripts`](https://mritchielab.github.io/FLAMES/reference/scmixology_lib10_transcripts.md)
  : scMixology long-read transcript counts - sample 2

- [`scmixology_lib90`](https://mritchielab.github.io/FLAMES/reference/scmixology_lib90.md)
  : scMixology short-read gene counts - sample 1

- [`find_bin()`](https://mritchielab.github.io/FLAMES/reference/find_bin.md)
  : Find path to a binary Wrapper for Sys.which to find path to a binary

## Pipelines (deprecated)

- [`sc_long_pipeline()`](https://mritchielab.github.io/FLAMES/reference/sc_long_pipeline.md)
  : Pipeline for Single Cell Data (deprecated)
- [`sc_long_multisample_pipeline()`](https://mritchielab.github.io/FLAMES/reference/sc_long_multisample_pipeline.md)
  : Pipeline for Multi-sample Single Cell Data (deprecated)
- [`bulk_long_pipeline()`](https://mritchielab.github.io/FLAMES/reference/bulk_long_pipeline.md)
  : Pipeline for bulk long read RNA-seq data processing (deprecated)
