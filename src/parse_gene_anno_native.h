#include <unordered_map>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <list>
#include <Rcpp.h>

#include "ParseGFF3.hpp"
#include "Pos.h"
#include "StartEndPair.hpp"

#ifndef GFF_DATA
#define GFF_DATA

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

        Rcpp::List
        to_R();

        void
        from_R(Rcpp::List list);

        std::unordered_map<std::string, std::vector<StartEndPair>>
        remove_transcript_duplicates(bool update_transcript_dict);
};

#endif

Rcpp::List
parse_gff_or_gtf_R(std::string filename);

GFFData
parse_gff_or_gtf(std::string filename);

GFFData
parse_gtf_tree(std::string filename);

GFFData
parse_gff_tree(std::string filename);