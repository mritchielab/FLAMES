#ifndef GROUP_BAM2ISOFORM_H
#define GROUP_BAM2ISOFORM_H

#include <vector>
#include <string>
#include <unordered_map>
#include <map>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <Rcpp.h>

#include "../classes/Config.h"
#include "../classes/GeneBlocks.h"
#include "../classes/Isoforms.h"
#include "../classes/StartEndPair.h"
#include "../classes/Pos.h"
#include "../classes/IsoKey.h"
#include "../classes/BamRecord.h"
#include "../utility/junctions.h"
// don't need to declare static function inside header
// static int
// fetch_function(const bam1_t *b, void *data);

void
bam_read (std::string bam_in, std::string chr, int s, int e);

std::vector<StartEndPair>
get_blocks(BAMRecord record);

void
minimal_group_bam2isoform (
    std::string bam_in, 
    std::string out_gff3, 
    std::string out_stat, 
    std::unordered_map<std::string, std::vector<GeneBlocks>>    * chr_to_blocks, 
    std::unordered_map<std::string, std::vector<StartEndPair>>  * gene_dict, 
    std::unordered_map<std::string, Junctions>                  * transcript_to_junctions,
    std::unordered_map<std::string, Pos>                        * transcript_dict,
    std::string fa_f,
    IsoformParameters isoform_parameters,
    std::string raw_gff3
);
void
group_bam2isoform (
    std::string bam_in, 
    std::string out_gff3, 
    std::string out_stat, 
    std::unordered_map<std::string, std::vector<GeneBlocks>>    * chr_to_blocks, 
    std::unordered_map<std::string, std::vector<StartEndPair>>  * gene_dict, 
    std::unordered_map<std::string, Junctions>                  * transcript_to_junctions,
    std::unordered_map<std::string, Pos>                        * transcript_dict,
    std::string fa_f,
    IsoformParameters isoform_parameters,
    std::string raw_gff3
);

#endif // GROUP_BAM2ISOFORM_H
