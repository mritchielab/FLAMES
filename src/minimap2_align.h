#include <string>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <Rcpp.h>

void
minimap2_align_cpp(std::string mm2_prog_path, std::string fa_file, std::string fq_in, std::string sam_out, bool no_flank, std::string bed12_junc);