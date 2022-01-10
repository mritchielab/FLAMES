#include <unordered_map>
#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include <iostream>

#include "../classes/BamRecord.h"

#include "../utility/bam.h"

#include "../main-functions/group_bam2isoform.h"

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

struct ReadDictEntry {
    std::string tr;
    float       AS_tag;
    float       tr_cov;
    float       length;
    int         quality;
};

#endif

std::unordered_map<std::string, int>
file_to_map(std::string filename);


void
read_entire_bam
(
    std::string bam_in, std::string log_out
);


int
query_len(std::string cigar_string, bool hard_clipping);


RealignedBamData
parse_realigned_bam
(
    std::string bam_in,
    std::string fa_idx_f,
    int         min_sup_reads,
    float       min_tr_coverage,
    float       min_read_coverage,
    std::unordered_map<std::string, std::string> kwargs
);

int
query_len(std::string cigar_string, bool hard_clipping=false);

std::unordered_map<std::string, std::string>
make_bc_dict(std::string bc_anno);