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

#ifndef ISOFORM_OBJECTS
#define ISOFORM_OBJECTS

struct IsoformObjects
{
    std::unordered_map<std::string, Pos>
    transcript_dict;
    
    std::unordered_map<std::string, Pos>
    transcript_dict_iso;
};

#endif

Rcpp::List
isoform_objects_to_R(IsoformObjects * isoform_objects);

IsoformObjects
isoform_objects_from_R(Rcpp::List list);

Rcpp::List
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