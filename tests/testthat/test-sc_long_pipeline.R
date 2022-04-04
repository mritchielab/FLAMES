test_that("sc_long_pipeline", {
  testthat::expect_s4_class(sc_long_pipeline(
    genome_fa = "/stornext/General/data/user_managed/grpu_mritchie_1/LuyiTian/Index/GRCm38.primary_assembly.genome.fa",
    fastq = "/stornext/General/data/user_managed/grpu_mritchie_1/Changqing/FLAMES_datasets/MuSC/sample/rps24.fastq",
    annot = "/stornext/General/data/user_managed/grpu_mritchie_1/LuyiTian/Index/gencode.vM24.annotation.gtf",
    outdir = "/home/users/allstaff/wang.ch/tests_FLAMES/sc_long_pipeline",
    match_barcode = FALSE,
    has_UMI = TRUE,
    minimap2_dir = "/home/users/allstaff/wang.ch/.conda/envs/FLAMES/bin",
    config_file = "/stornext/General/data/user_managed/grpu_mritchie_1/Changqing/config_sclr_nanopore_5end.json",
    seed = 2022
  ), "SingleCellExperiment")
})

test_that("sc_long_pipeline_bambu", {
  testthat::expect_s4_class(sc_long_pipeline(
    genome_fa = "/stornext/General/data/user_managed/grpu_mritchie_1/LuyiTian/Index/GRCm38.primary_assembly.genome.fa",
    fastq = "/stornext/General/data/user_managed/grpu_mritchie_1/Changqing/FLAMES_datasets/MuSC/sample/rps24.fastq",
    annot = "/stornext/General/data/user_managed/grpu_mritchie_1/LuyiTian/Index/gencode.vM24.annotation.gtf",
    outdir = "/home/users/allstaff/wang.ch/tests_FLAMES/sc_long_pipeline_bambu",
    match_barcode = FALSE,
    has_UMI = TRUE,
    isoform_id_bambu = TRUE,
    in_bam = "/stornext/General/data/user_managed/grpu_mritchie_1/Changqing/FLAMES_datasets/MuSC/sample/rps24.bam",
    minimap2_dir = "/home/users/allstaff/wang.ch/.conda/envs/FLAMES/bin",
    config_file = "/stornext/General/data/user_managed/grpu_mritchie_1/Changqing/config_sclr_nanopore_5end.json",
    seed = 2022
  ), "SingleCellExperiment")
})