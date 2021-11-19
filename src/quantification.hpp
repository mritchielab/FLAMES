#include "config.h"
#include "parse_realigned_bam.hpp"
#include "write_tr_to_csv.hpp"
#include "annotate_filter_gff.hpp"
#include "find_isoform.h"

void
quantification_cpp
(
    Rcpp::List  config_list,
    std::string realign_bam,
    std::string transcript_fa_idx,
    Rcpp::List  isoform_objects_list,
    std::string tr_cnt_csv,
    std::string tr_badcov_cnt_csv,
    std::string isoform_gff3,
    std::string annot,
    std::string isoform_gff3_f,
    std::string FSM_anno_out
);