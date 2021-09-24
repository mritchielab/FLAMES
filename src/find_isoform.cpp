#include "find_isoform.h"

void
find_isoform
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
    
    std::unordered_map<std::string, std::vector<std::string>>
    chr_to_gene;
    std::unordered_map<std::string, Pos>
    transcript_dict;
    std::unordered_map<std::string, std::vector<std::string>>
    gene_to_transcript;
    std::unordered_map<std::string, std::vector<std::pair<int, int>>>
    transcript_to_exon;

    
}