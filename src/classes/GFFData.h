#ifndef GFF_DATA
#define GFF_DATA

#include <unordered_map>
#include <string>
#include <vector>

#include "ParseGFF3.hpp"
#include "Pos.h"
#include "StartEndPair.hpp"

class GFFData
{
    public:
        std::unordered_map<std::string, std::vector<std::string>>
        chr_to_gene;
    
        std::unordered_map<std::string, Pos>
        transcript_dict;

        std::unordered_map<std::string, std::vector<std::string>>
        gene_to_transcript;

        std::unordered_map<std::string, std::vector<StartEndPair>>
        transcript_to_exon;

        Rcpp::List to_R();

        void from_R(Rcpp::List list);

        void remove_transcript_duplicates(bool update_transcript_dict);
};