#include "minimap2_align.h"

// [[Rcpp::export]]
void
minimap2_align_cpp
(
    std::string mm2_prog_path,
    std::string fa_file,
    std::string fq_in,
    std::string sam_out,
    bool no_flank,
    std::string bed12_junc
)
{
    /*
        calls minimap2 from the command line to align to genome, given all the necessary parameters
        default command is:
        {mm2_prog_path}/minimap2 -ax splice -t 12  -k14 --secondary=no {fa_file} {fq_in} -o {sam_out}
    */

    std::string junc_cmd;
    if (bed12_junc != "") {
        junc_cmd = "--junc-bed " + bed12_junc + " --junc-bonus 1";
    } else {
        junc_cmd = "";
    }

    std::string flank_cmd;
    if (no_flank) {
        flank_cmd = "--splice-flank=no";
    } else {
        flank_cmd = "";
    }

    // form the command from all the separate parts
    std::stringstream align_cmd;
    align_cmd << mm2_prog_path << "minimap2" << " -ax splice -t 12 " 
            << junc_cmd << " " 
            << flank_cmd << " -k14 --secondary=no " 
            << fa_file << " " 
            << fq_in << " -o " << sam_out;
    // call it
    system(align_cmd.str().c_str());
}

// [[Rcpp::export]]
void
minimap2_tr_align_cpp
(
    std::string mm2_prog_path,
    std::string fa_file,
    std::string fq_in,
    std::string sam_out
)
{
    /* calls minimap2 to align to transcript */

    std::stringstream align_cmd;
    align_cmd << mm2_prog_path << "minimap2" 
            << " -ax map-ont -p 0.9 --end-bonus 10 -N 3 -t 12 " 
            << fa_file << " " 
            << fq_in << " -o " 
            << sam_out;
    // call it
    system(align_cmd.str().c_str());
}