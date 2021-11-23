#ifndef REFERENCE_DICT
#define REFERENCE_DICT

#include <string>
#include <vector>
#include <unordered_map>
#include "Pos.h"
#include "StartEndPair.hpp"

// a quick struct for a reference dictionary
struct ReferenceDict {
    std::unordered_map<std::string, std::vector<std::string>>
    chr_to_gene;
    std::unordered_map<std::string, Pos>
    transcript_dict;
    std::unordered_map<std::string, std::vector<std::string>>
    gene_to_transcript;
    std::unordered_map<std::string, std::vector<StartEndPair>>
    transcript_to_exon;
};

#endif