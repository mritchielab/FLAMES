#include "find_isoform.h"

// [[Rcpp::export]]
void
find_isoform_cpp
(
    std::string gff3, 
    std::string genome_bam, 
    std::string isoform_gff3, 
    std::string tss_test_stat, 
    std::string genomefa, 
    std::string transcript_fa, 
    int downsample_ratio, 
    Rcpp::List config_list, 
    std::string raw_splice_isoform
)
{
    std::cout << "### Read Gene Annotations";

    // first, extract from the gff file
    Rcpp::List gene_anno = parse_gff_tree(gff3.c_str());
    
    // convert the chr_to_gene Rcpp::List to a std::unordered_map
    std::unordered_map<std::string, std::vector<std::string>>
    chr_to_gene;
    // for (const auto & [key, val] : gene_anno["chr_to_gene"]) {
    //     chr_to_gene[Rcpp::as(key)] = Rcpp::as(val);
    // }
    
    std::unordered_map<std::string, Pos>
    transcript_dict;

    // for (const auto & [key, val] : gene_anno["chr_to_gene"]) {
    //     chr_to_gene[Rcpp::as(key)] = val;
    // }
    
    std::unordered_map<std::string, std::vector<std::string>>
    gene_to_transcript;
    std::unordered_map<std::string, std::vector<std::pair<int, int>>>
    transcript_to_exon;

    

    
}