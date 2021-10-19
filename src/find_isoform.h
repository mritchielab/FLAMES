#include <string>
#include <iostream>
#include <unordered_map>
#include <vector>
#include <Rcpp.h>
#include "config.h"
#include "misc.h"
#include "parse_gene_anno_native.h"
#include "junctions.h"
#include "gff3_to_fa.hpp"
#include "group_bam2isoform.h"
#include "ReferenceDict.hpp"

void
find_isoform_cpp
(
    std::string gff3, 
    std::string genome_bam, 
    std::string isoform_gff3, 
    std::string tss_test_stat, 
    std::string genomefa, 
    std::string transcript_fa, 
    int         downsample_ratio, 
    Rcpp::List  config_list, 
    std::string raw_splice_isoform
);