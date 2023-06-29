#ifndef GFFDATA_H
#define GFFDATA_H

#include <unordered_map>
#include <vector>
#include <string>
#include <Rcpp.h>

#include "Pos.h"            // Pos
#include "types.h"          // transcriptvector, exon

struct GFFData
{
    std::unordered_map<std::string, std::vector<std::string>>
    chr_to_gene;

    std::unordered_map<std::string, Pos>
    transcript_dict;

    std::unordered_map<std::string, transcriptvector>
    gene_to_transcript;

    std::unordered_map<std::string, std::vector<exon>>
    transcript_to_exon;

    // Rcpp::List to_R();

    // void
    // from_R(Rcpp::List list);

    // void
    // log(std::string filename);

    bool is_empty() const {
        return (chr_to_gene.size() + transcript_dict.size() + gene_to_transcript.size() + transcript_to_exon.size()) > 0;
    }
};

#endif