#ifndef GFFDATA_H
#define GFFDATA_H

#include <unordered_map>
#include <vector>
#include <string>
#include <Rcpp.h>
#include <fstream>

#include "Pos.h"
#include "StartEndPair.h"

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

        void
        from_R(Rcpp::List list);

        void
        removeTranscriptDuplicates(bool updateTranscriptDict=true);

        void
        log(std::string filename);
};

#endif