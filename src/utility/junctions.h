/*
  Anything related to junctions are goes in here
*/
#ifndef JUNCTIONS
#define JUNCTIONS

#include <vector>
#include <string>
#include <set>
#include <unordered_map>

#include "GeneBlocks.h"
#include "StartEndPair.hpp"

struct DoubleJunctions {
  /*
    a struct used for holding the junctions information
  */
  std::vector<int> left;
//   std::vector<int> junctions;
  std::vector<int> right;
};

//struct SingleJunction {
struct Junctions {
  /*
    a struct for junctions with just one value in left and right
  */
  int left;
  std::vector<int> junctions;
  int right;
};

int 
take_closest (std::vector<int> list, int num);

Junctions 
blocks_to_junctions (std::vector<StartEndPair> blocks);

DoubleJunctions 
get_TSS_TES_site
(
    std::unordered_map<std::string, Junctions> * transcript_to_junctions,
    const std::vector<std::string> * tr_list
);

std::set<int> 
get_splice_site (std::unordered_map<std::string, Junctions> transcript_to_junctions, std::vector<std::string> tr_list);

bool
is_exon_similar
(
    const std::vector<StartEndPair> &exon1, 
    const std::vector<StartEndPair> &exon2, 
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
    std::unordered_map<std::string, std::vector<std::string>>   &gene_to_transcript,
    const std::unordered_map<std::string, std::vector<StartEndPair>>  &transcript_to_exon,
    int threshold
);

std::unordered_map<std::string, std::vector<GeneBlocks>>
get_gene_blocks
(
    std::unordered_map<std::string, std::vector<StartEndPair>>  * gene_dict,
    std::unordered_map<std::string, std::vector<std::string>>   * chr_to_gene,
    std::unordered_map<std::string, std::vector<std::string>>   * gene_to_transcript
);

#endif
