#include <string>
#include <vector>
#include <unordered_map>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <limits>
#include <set>

// #include "parse_gene_anno_native.h"
#include "classes/GeneAnnoParser/GeneAnnoParser.h"

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