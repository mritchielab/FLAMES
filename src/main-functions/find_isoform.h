#ifndef FIND_ISOFORMS_H
#define FIND_ISOFORMS_H

#include <unordered_map>
#include <string>
#include <vector>
#include <Rcpp.h>

#include "../classes/Config.h"
#include "../classes/GFFData.h"
#include "../classes/ReferenceDict.h"
#include "../classes/GeneAnnoParser/GeneAnnoParser.h"
#include "../classes/Pos.h"

#include "../utility/misc.h"
#include "../utility/junctions.h"

#include "../file-handling/gff3_to_fa.h"

#include "group_bam2isoform.h"


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


#endif // FIND_ISOFORMS_H