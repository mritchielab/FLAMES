#ifndef JUNCTIONS_H
#define JUNCTIONS_H

#include <vector>
#include <string>
#include <unordered_map>
#include <set>

#include "StartEndPair.h"
#include "types.h"
#include "GeneBlocks.h"

struct DoubleJunctions {
    /*
        a struct used for holding the junctions information
    */
    std::vector<int> left;
    // std::vector<int> junctions;
    std::vector<int> right;
};

//struct SingleJunction {
struct Junctions {
    /*
        a struct for junctions with just one value in left and right
    */
    int left;
    std::vector<int> junctions;
    int right;
};

struct InvJunctions {
    /*
        a struct for junctions with just one value in left and right
    */
    int right;
    std::vector<int> junctions;
    int left;
};

Junctions 
blocks_to_junctions (const std::vector<StartEndPair> &blocks);

std::unordered_map<std::string, Junctions>
map_transcripts_to_junctions(
    const std::unordered_map<std::string, std::vector<exon>> &transcript_to_exons
);

std::unordered_map<std::string, transcriptvector>
remove_similar_tr(
    const std::unordered_map<std::string, transcriptvector> &gene_to_transcript,
    const std::unordered_map<std::string, std::vector<exon>>  &transcript_to_exons,
    int threshold
);

bool
is_exon_similar (
    const std::vector<exon> &exon1, 
    const std::vector<exon> &exon2, 
    int threshold);

std::unordered_map<std::string, std::vector<exon>>
get_gene_flat(
    const std::unordered_map<std::string, transcriptvector> &gene_to_transcript,
    const std::unordered_map<std::string, std::vector<exon>> &transcript_to_exons);

std::unordered_map<std::string, std::vector<GeneBlocks>>
get_gene_blocks(
    const std::unordered_map<std::string, std::vector<exon>> &gene_dict, 
    const std::unordered_map<std::string, std::vector<std::string>> &chr_to_gene,
    const std::unordered_map<std::string, transcriptvector> &gene_to_transcript);

DoubleJunctions 
get_TSS_TES_site(
    const std::unordered_map<std::string, Junctions> &transcript_to_junctions,
    const std::vector<std::string> &tr_list
);

int 
take_closest (const std::vector<int> &list, int num);

std::set<int> 
get_splice_site (
    const std::unordered_map<std::string, Junctions> &transcript_to_junctions, 
    const std::vector<std::string> &tr_list
);

std::vector<int>
find_best_splice_chain(const std::vector<int> &raw_iso, const std::vector<std::vector<int>> &junction_list, int MAX_DIST);

int
exon_overlap(const std::vector<exon> &exons1, const std::vector<exon> &exons2);
int
exon_overlap (const std::vector<int> &exons1, const std::vector<exon> &exons2);

std::vector<exon>
pairwiseStartEndPair (const std::vector<int> &input);

bool
check_exon_subset(const std::vector<int> &s1, const std::vector<int> &s2, int MAX_TOLERANCE);

float
get_exon_sim_pct(const std::vector<int> &exons1, const std::vector<int> &exons2);
#endif // JUNCTIONS_H