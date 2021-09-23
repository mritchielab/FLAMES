#include <vector>
#include <string>
#include <map>

#ifndef GENEBLOCKS
#define GENEBLOCKS

class GeneBlocks
{
  /*
    a class for storing info about individual gene blocks
    can be updated with add_gene
  */
  public:
    GeneBlocks(int start, int end, std::vector<std::string> transcript_list, std::string a_gene);
    
    int start, end;
    std::vector<std::string> 
    transcript_list;
    std::map<std::string, std::vector<std::string>>
    gene_to_transcript;
    std::map<std::string, std::vector<std::string>>
    gene_to_tr;

    void add_gene(int start, int end, std::vector<std::string> transcript_list, std::string a_gene);
};

#endif