#ifndef ANNOTATE_FILTER_GFF_H
#define ANNOTATE_FILTER_GFF_H

#include <string>
#include <unordered_map>
#include <vector>

#include "../classes/Pos.h"
#include "../classes/StartEndPair.h"
void
annotate_filter_gff
(
    std::string isoform_gff,
    std::string ref_gff,
    std::string isoform_out,
    std::string anno_out,
    std::unordered_map<std::string, int> tr_count,
    int min_sup_reads,
    bool verbose=true
);

void
annotate_full_splice_match
(
    std::unordered_map<std::string, std::vector<StartEndPair>> transcript_to_exon,
    std::unordered_map<std::string, std::vector<StartEndPair>> transcript_to_exon_ref,
    std::unordered_map<std::string, Pos> transcript_dict,
    std::unordered_map<std::string, Pos> transcript_dict_ref,
	std::string anno_out,
    std::unordered_map<std::string, int> tr_count,
    int min_sup_reads
);

#endif // ANNOTATE_FILTER_GFF_H