#ifndef GENEBLOCKS_H
#define GENEBLOCKS_H

#include <vector>
#include <string>
#include <map>
#include <algorithm>

#include "types.h"

struct GeneBlocks
{
    int start, end;
    transcriptvector transcript_list;
    std::map<std::string, transcriptvector> gene_to_transcript;

    GeneBlocks(int _start, int _end, const transcriptvector &_transcript_list, const std::string &a_gene)
        : start{_start}, end{_end}, transcript_list{_transcript_list}, gene_to_transcript{{a_gene, _transcript_list}} 
        {}

    void add_gene(int _start, int _end, const transcriptvector &_transcript_list, const std::string &a_gene)
    {
        end = std::max(end, _end);
        transcript_list.insert(transcript_list.end(), _transcript_list.begin(), _transcript_list.end());
        gene_to_transcript[a_gene] = transcriptvector(transcript_list);
    }
};
#endif