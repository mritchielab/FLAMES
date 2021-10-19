#include <string>
#include <map>
#include <vector>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <Rcpp.h>

#include "GeneBlocks.h"
#include "junctions.h"
#include "misc.h"
#include "Isoforms.h"
#include "bam.h"

using namespace Rcpp;


static int
fetch_function(const bam1_t *b, void *data);

void
bam_read (std::string bam_in, int s, int e);

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
    Config config,
    std::string raw_gff3
);
