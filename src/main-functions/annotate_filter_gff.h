#ifndef ANNOTATE_FILTER_GFF_H
#define ANNOTATE_FILTER_GFF_H

#include <string>
#include <unordered_map>

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

#endif // ANNOTATE_FILTER_GFF_H