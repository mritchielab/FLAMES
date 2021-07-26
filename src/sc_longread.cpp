#include <Rcpp.h>
#include <R.h>
#include <iostream>
#include <sstream>
#include <fstream>

#include <unordered_map>

using namespace Rcpp;

//' Get Gene Blocks
//'
//' @description
//'
//' @param gene_dict std::unordered_map of ?
//' @param chr_to_gene std::unordered_map of  - output of parse_gff_tree
//' @param gene_to_transcript std::unordered_map of - output of parse_gff_tree
// void get_gene_blocks(std::unordered_map<key, value> gene_dict, 
//                     std::unordered_map<std::string, StringVector> chr_to_gene, 
//                     std::unordered_map<key, value> gene_to_transcript) {
//     // output map of chromosomes to gene blocks
//     std::unordered_map<key, value> chr_to_blocks;
// }
void get_gene_blocks(?? gene_dict, 
                    std::unordered_map<String, std::unordered_map<String, bool>> chr_to_gene, 
                    std::unordered_map<String, std::unordered_map<String, bool>> gene_to_transcript) {
                        
                    }