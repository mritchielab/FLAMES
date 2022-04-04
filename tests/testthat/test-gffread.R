test_that("gffread", {
  test_out <- "/home/users/allstaff/wang.ch/tests_FLAMES/gffread_cpp.fa"
  compare_file <- "/stornext/General/data/user_managed/grpu_mritchie_1/Changqing/FLAMES_datasets/MuSC/sample/gffread.fa"
  gffread_cpp(
    "/stornext/General/data/user_managed/grpu_mritchie_1/LuyiTian/Index/GRCm38.primary_assembly.genome.fa",
    test_out,
    "/stornext/General/data/user_managed/grpu_mritchie_1/Changqing/FLAMES_datasets/MuSC/sample/gffread.gff3"
  )
  testthat::expect_true(all.equal(readLines(test_out), readLines(compare_file)))
})