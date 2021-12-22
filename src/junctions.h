/*
  Anything related to junctions are goes in here
*/

#include <stdio.h>
#include <iostream>
#include <map>
#include <vector>
#include <any>
#include <set>
#include <list>

#include "GeneBlocks.h"
#include "ParseGFF3.hpp"
#include "StartEndPair.hpp"

#ifndef JUNCTIONS
#define JUNCTIONS

struct Junctions {
  /*
    a struct used for holding the junctions information
  */
  std::vector<int> left;
  std::vector<int> junctions;
  std::vector<int> right;
};

struct SingleJunction {
  /*
    a struct for junctions with just one value in left and right
  */
  int left;
  std::vector<int> junctions;
  int right;
};

#endif

int 
take_closest (std::vector<int> list, int num);

Junctions 
blocks_to_junctions (std::vector<StartEndPair> blocks);

Junctions 
get_TSS_TES_site
(
    std::unordered_map<std::string, Junctions> * transcript_to_junctions,
    const std::vector<std::string> * tr_list
);

std::set<int> 
get_splice_site (std::unordered_map<std::string, Junctions> transcript_to_junctions, std::vector<std::string> tr_list);

int
is_exon_similar
(
    std::vector<StartEndPair> * exon1, 
    std::vector<StartEndPair> * exon2, 
    int threshold
);

std::unordered_map<std::string, std::vector<StartEndPair>>
get_gene_flat(
    std::unordered_map<std::string, std::vector<std::string>>   * gene_to_transcript,
    std::unordered_map<std::string, std::vector<StartEndPair>>  * transcript_to_exon
);

void
remove_similar_tr
(
    std::unordered_map<std::string, std::vector<std::string>>   * gene_to_transcript,
    std::unordered_map<std::string, std::vector<StartEndPair>>  * transcript_to_exon,
    int threshold
);

std::unordered_map<std::string, std::vector<GeneBlocks>>
get_gene_blocks
(
    std::unordered_map<std::string, std::vector<StartEndPair>>  * gene_dict,
    std::unordered_map<std::string, std::vector<std::string>>   * chr_to_gene,
    std::unordered_map<std::string, std::vector<std::string>>   * gene_to_transcript
);

void
junctions_print(Junctions junctions);