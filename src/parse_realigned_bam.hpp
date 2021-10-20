#include <unordered_map>
#include <string>
#include <fstream>
#include <vector>
#include <sstream>

void
parse_realigned_bam(
    std::string bam_in,
    std::string fa_idx_f,
    std::string min_sup_reads,
    std::string min_tr_coverage,
    std::string min_read_coverage,
    std::string kwargs
);

