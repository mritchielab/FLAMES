#ifndef MATCH_CELL_BARCODE_H
#define MATCH_CELL_BARCODE_H

#include <vector>
#include <utility>
#include <unordered_map>
#include <string>

#include <Rcpp.h>

#include "../utility/edit_dist.h"
#include "../utility/ssw/ssw_cpp.h"
#include "../utility/fastq_utils.h"

std::string
join_path
(
    const std::string p1, const std::string p2
);

int
find_polyT
(
    std::string &seq, 
    int start_pos
);

char complement(char n);

void
rc
(
    std::string &seq
);

int
getdir
(
    const char dir[], 
    std::vector<std::string> &files
);

void
get_bc_anno
(
    std::string fn, 
    std::unordered_map<std::string, std::string> &barcode_dict, 
    std::vector<std::string> &barcode_list
);

bool
sortbysec_dec
(
    const std::pair<int, int> &a,
    const std::pair<int, int> &b
);

std::pair<int, int> 
get_bc_range
(
    std::string fqn, 
    int max_reads, 
    const std::string lf_seq, 
    int MAX_DIST
);

int
get_clo_idx
(
    const char *seq_ptr, 
    int64_t *al, 
    std::vector<int64_t *> &bc_list_ptr, 
    int max_dist
);

int
get_hm_idx
(
    std::string &q_seq, 
    std::vector<std::string> &barcode_list, 
    int max_dist
);

Rcpp::List
match_cell_barcode
(
    Rcpp::String fastq_dir, 
    Rcpp::String stats_file, 
    Rcpp::String out_fastq, 
    Rcpp::String ref_csv, 
    int MAX_DIST, 
    int UMI_LEN
);

#endif // MATCH_CELL_BARCODE_H
