#ifndef FIND_ISOFORMS_H
#define FIND_ISOFORMS_H

#include <unordered_map>
#include <string>
#include <vector>
#include <Rcpp.h>

#include "../classes/Config.h"
#include "../classes/GFFData.h"
#include "../utility/junctions.h"
#include "../file-handling/gff3_to_fa.h"
#include "../main-functions/group_bam2isoform.h"
#include "../classes/ReferenceDict.h"
#include "../classes/Pos.h"

struct IsoformObjects
{
    std::unordered_map<std::string, Pos>
    transcript_dict;
    
    std::unordered_map<std::string, Pos>
    transcript_dict_iso;
};


Rcpp::List
isoform_objects_to_R(IsoformObjects * isoform_objects);

IsoformObjects
isoform_objects_from_R(Rcpp::List list);

Rcpp::List
find_isoform
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

void
log_params
(
    std::unordered_map<std::string, std::vector<GeneBlocks>> *
    chr_to_blocks,

    std::unordered_map<std::string, std::vector<StartEndPair>> *
    gene_dict,

    std::unordered_map<std::string, Junctions> *
    transcript_to_junctions,

    std::unordered_map<std::string, Pos> *
    transcript_dict,

    std::string
    filename="group_bam2isoform_param_cpp.txt"
);

#endif // FIND_ISOFORMS_H