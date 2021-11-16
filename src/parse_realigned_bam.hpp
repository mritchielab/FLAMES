#include <unordered_map>
#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include <iostream>

#include "bam.h"
#include "group_bam2isoform.h"

#ifndef REALIGNED_BAM_DATA
#define REALIGNED_BAM_DATA

/*  a struct to package up the output of parse_realigned_bam
*/
struct RealignedBamData {
    std::unordered_map<std::string, std::unordered_map<std::string, std::vector<std::string>>>
    bc_tr_count_dict;

    std::unordered_map<std::string, std::unordered_map<std::string, std::vector<std::string>>>
    bc_tr_badcov_count_dict;

    std::vector<std::string>
    tr_kept;
};

#endif

RealignedBamData
parse_realigned_bam(
    std::string bam_in,
    std::string fa_idx_f,
    std::string min_sup_reads,
    std::string min_tr_coverage,
    std::string min_read_coverage,
    std::string kwargs
);


int
query_len(std::string cigar_string, bool hard_clipping=false);

std::unordered_map<std::string, std::string>
make_bc_dict(std::string bc_anno);