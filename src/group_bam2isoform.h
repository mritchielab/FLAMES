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
#include "StartEndPair.hpp"
#include "cigars.h"

using namespace Rcpp;

#ifndef REC
#define REC

/*  a struct for handling each record in a bam file
*/
struct Record {
    std::vector<CigarPair>
    cigar;
    std::string
    cigar_string;

    int
    reference_start;
    int
    reference_end;
    bool
    is_reverse;
};

#define BAM_CMATCH      0   // CIGAR character for matching
#define BAM_CDEL        2   
#define BAM_CREF_SKIP   3
#define BAM_CEQUAL      7
#define BAM_CDIFF       8

#endif

static int
fetch_function(const bam1_t *b, void *data);

void
bam_read (std::string bam_in, int s, int e);

std::vector<StartEndPair>
get_blocks(Record record);

void
group_bam2isoform (
    std::string bam_in, 
    std::string out_gff3, 
    std::string out_stat, 
    std::unordered_map<std::string, std::vector<GeneBlocks>>    * chr_to_blocks, 
    std::unordered_map<std::string, std::vector<StartEndPair>   * gene_dict, 
    std::unordered_map<std::string, Junctions>                  * transcript_to_junctions,
    std::unordered_map<std::string, Pos>                        * transcript_dict,
    std::string fa_f,
    IsoformParameters isoform_parameters,
    std::string raw_gff3 = ""
);

