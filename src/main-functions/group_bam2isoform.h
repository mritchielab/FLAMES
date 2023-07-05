#ifndef GROUP_BAM2ISOFORM_H
#define GROUP_BAM2ISOFORM_H

#include <string>
#include <vector>
#include <unordered_map>

#include "../classes/GFFData.h"
#include "../classes/StartEndPair.h"
#include "../classes/junctions.h"
#include "../classes/Pos.h"
#include "../classes/GeneBlocks.h"

void group_bam2isoform(
    const std::string &bam_in,
    const std::string &out_gff3,
    const std::string &out_stat,
    const std::unordered_map<std::string, std::vector<StartEndPair>>  	&gene_dict, 
    const std::unordered_map<std::string, Junctions>           	 		&transcript_to_junctions,
    const std::unordered_map<std::string, Pos>                        	&transcript_dict,
    const std::unordered_map<std::string, std::vector<GeneBlocks>>      &chr_to_blocks,
    const std::string &fa_file,
    const Rcpp::List &isoform_parameters,
    const std::string &raw_gff3
);

#endif // GROUP_BAM2ISOFORM_H