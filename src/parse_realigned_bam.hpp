#include <unordered_map>
#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include <iostream>

#include "bam.h"
#include "group_bam2isoform.h"

void
parse_realigned_bam(
    std::string bam_in,
    std::string fa_idx_f,
    std::string min_sup_reads,
    std::string min_tr_coverage,
    std::string min_read_coverage,
    std::string kwargs
);


std::unordered_map<std::string, std::string>
make_bc_dict(std::string bc_anno);