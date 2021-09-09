#include <string>
#include <map>
#include <vector>
#include <fstream>

#include "GeneBlocks.h"
#include "junctions.h"
#include "misc.h"
#include "Isoforms.h"
#include "bam.h"

void
group_bam2isoform (
    std::string bam_in, 
    std::string out_gff3, 
    std::string out_stat, 
    std::map<std::string, 
    std::vector<GeneBlocks>> chr_to_blocks, 
    std::map<std::string, std::vector<std::pair<int, int>>> gene_dict, 
    std::map<std::string, Junctions> transcript_to_junctions,
    std::map<std::string, Pos> transcript_dict,
    std::string fa_f,
    std::map<std::string, int> config,
    std::string raw_gff3
);
