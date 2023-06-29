#include "junctions.h"

#include <vector>
#include <unordered_map>
#include <string>
#include <unordered_set>
#include <utility>
#include <set>

#include "StartEndPair.h"
#include "types.h"
#include "../utility/utility.h"
#include "GeneBlocks.h"

/*
    Converts between a list of StartEndPair objects (blocks)
    into a Junctions object, representing a list of all positions,
    bounded by start of the first block, and end of the last block
*/
Junctions 
blocks_to_junctions (const std::vector<StartEndPair> &blocks)
{
    Junctions output;

    if ((int)blocks.size() < 1) return output;

	output.left = blocks.front().start;
	output.right = blocks.back().end;

    if ((int)blocks.size() == 1) return output;

    // Fill in the Junctions.junctions vector
    for (int i = 1; i < (int)blocks.size(); i++) {
        output.junctions.push_back(blocks[i - 1].end);
        output.junctions.push_back(blocks[i].start);
    }

    return output; 
}

/*
    Convert all exons for a transcript into a Junctions object,
    and associate with that transcript in a new unordered map
*/
std::unordered_map<std::string, Junctions>
map_transcripts_to_junctions(
    const std::unordered_map<std::string, std::vector<exon>> &transcript_to_exons)
{
    std::unordered_map<std::string, Junctions> transcript_to_junctions;
    for (const auto &[tr, ex] : transcript_to_exons) {
        transcript_to_junctions[tr] = blocks_to_junctions(ex);
    }

    return transcript_to_junctions;
}

/* 
    takes two vectors of exons,
    checks that they are similar enough within a given threshold
*/
bool
is_exon_similar (
    const std::vector<exon> &exon1, 
    const std::vector<exon> &exon2, 
    int threshold)
{
    // make sure they're the same size
    if (exon1.size() != exon2.size()) return false;

    // then add up the difference
    int diff = 0;
    for (int i = 0; i < (int)exon1.size(); ++i) {
        diff += std::abs(exon1.at(i).start - exon2.at(i).start) 
                + std::abs(exon1.at(i).end - exon2.at(i).end);

        if (diff > threshold) return false;
    }

    return true;
}

/*
    Remove similar transcripts from gene_to_transcript based on
    using the resultant exon from transcript_to_exon as a comparator.
    Returns a new std::unordered_map<std::string, std::vector<std::string>> containing
    only the transcripts that aren't similar for each gene
*/
std::unordered_map<std::string, transcriptvector>
remove_similar_tr(
    const std::unordered_map<std::string, transcriptvector>   &gene_to_transcript,
    const std::unordered_map<std::string, std::vector<exon>>  &transcript_to_exon,
    int threshold)
{
    int dup_count = 0;
    int transcriptIdx1, transcriptIdx2;

    std::unordered_map<std::string, transcriptvector> out_gene_to_transcript;
    
    for (const auto &[gene, transcript] : gene_to_transcript) {
        // ignore anything short
        if (transcript.size() < 2) continue;

        // keep a set of all the duplicates
        std::unordered_set<int> dup_set;

        // iterate through every pair of transcripts t1 and t2 in the gene
        for (transcriptIdx1 = 0; transcriptIdx1 < (int)transcript.size() - 1; ++transcriptIdx1) {
            // Don't check for duplicates if already removing this transcript
            if (dup_set.count(transcriptIdx1)) continue;
            
            for (transcriptIdx2 = transcriptIdx1 + 1; transcriptIdx2 < (int)transcript.size(); ++transcriptIdx2) {
                // Don't check against duplicates if already removing this transcript
                if (dup_set.count(transcriptIdx2)) continue;
                
                if (is_exon_similar(
                    transcript_to_exon.at(transcript[transcriptIdx1]), 
                    transcript_to_exon.at(transcript[transcriptIdx2]), 
                    threshold)) 
                {
                    dup_set.insert(transcriptIdx2);
                    dup_count++;
                }
            }
        }
        
        // Create a new transcriptvector of the transcripts that weren't removed
        out_gene_to_transcript[gene] = {};
        for (int i = 0; i < (int)transcript.size(); i++) {
            // Add this transcript to the transcriptvector if it wasn't marked as a duplicate
            if (dup_set.count(i) == 0) {
                out_gene_to_transcript.at(gene).push_back(transcript.at(i));
            }
        }
    }

    return out_gene_to_transcript;
}

std::unordered_map<std::string, std::vector<exon>>
get_gene_flat(
    const std::unordered_map<std::string, transcriptvector> &gene_to_transcript,
    const std::unordered_map<std::string, std::vector<exon>> &transcript_to_exons)
{
    std::unordered_map<std::string, std::vector<exon>> gene_dict;

    for (const auto &[gene, transcripts] : gene_to_transcript) {
        std::vector<exon> exons;
        // Create a contiguous list of all the exons in transcript_to_exon[transcript]
        for (const std::string &transcript : transcripts) {
            const auto &insertExon = transcript_to_exons.at(transcript);
            exons.insert(exons.end(), insertExon.begin(), insertExon.end());
        }
        // just store the exon if there's only one
        if (exons.size() == 1) {
            gene_dict[gene] = {exon(exons[0])};
            continue;
        }

        // Organise the exons and merge
        ranges::sort(exons, StartEndPairCompare);
        std::vector<exon> merged_exons = {exons.at(0)};
        for (int i = 1; i < exons.size(); i++) {
            const exon &higher = exons.at(i);
            const exon &lower = merged_exons.back();

            if (higher.start <= lower.end) {
                int end = std::max(lower.end, higher.end);
                merged_exons.back().start = lower.start;
                merged_exons.back().end = end;
            } else {
                merged_exons.push_back(higher);
            }
        }
        gene_dict[gene] = merged_exons;
    }

    return gene_dict;
}

std::unordered_map<std::string, std::vector<GeneBlocks>>
get_gene_blocks(
    const std::unordered_map<std::string, std::vector<exon>> &gene_dict, 
    const std::unordered_map<std::string, std::vector<std::string>> &chr_to_gene,
    const std::unordered_map<std::string, transcriptvector> &gene_to_transcript)
{
    std::unordered_map<std::string, std::vector<GeneBlocks>> chr_to_blocks;

    for (const auto &[chr, genes] : chr_to_gene) {
        std::vector<std::pair<std::string, StartEndPair>> gene_list;
        for (const auto &gene : genes) {
            gene_list.push_back(
                {gene, 
                    StartEndPair(
                        gene_dict.at(gene).front().start, 
                        gene_dict.at(gene).back().end)
                }
            );
        }

        ranges::sort<std::pair<std::string, StartEndPair>>(gene_list, [](const auto &a, const auto &b) { return StartEndPairCompare(a.second, b.second); });

        auto first_gene = gene_list.front();
        chr_to_blocks[chr].push_back(GeneBlocks(
            first_gene.second.start,
            first_gene.second.end,
            gene_to_transcript.at(first_gene.first),
            first_gene.first
        ));

        if (gene_list.size() <= 1) {
            continue;
        }

        // Process the rest of the genes
        for (auto entry = gene_list.begin() + 1;  entry != gene_list.end(); entry++) {
            if (chr_to_blocks[chr].back().end > entry->second.start) {
                chr_to_blocks[chr].back().add_gene(
                    entry->second.start, 
                    entry->second.end, 
                    gene_to_transcript.at(entry->first),
                    entry->first
                );
            } else {
                chr_to_blocks[chr].push_back(GeneBlocks(
                    entry->second.start,
                    entry->second.end,
                    gene_to_transcript.at(entry->first),
                    entry->first
                ));
            }
        }
    }

    return chr_to_blocks;
}


DoubleJunctions 
get_TSS_TES_site(
    const std::unordered_map<std::string, Junctions> &transcript_to_junctions,
    const std::vector<std::string> &tr_list)
{
    DoubleJunctions all_site;

    for (const auto & t : tr_list) {
        if (all_site.left.size() > 0) {
            if (abs(take_closest(all_site.left, transcript_to_junctions.at(t).left) - transcript_to_junctions.at(t).left) > 5) {
                all_site.left.push_back(transcript_to_junctions.at(t).left);
            }
        } else {
            all_site.left.push_back(transcript_to_junctions.at(t).left);
        }

        if (all_site.right.size() > 0) {
            if (abs(take_closest(all_site.right, transcript_to_junctions.at(t).right) - transcript_to_junctions.at(t).right) > 5) {
                all_site.right.push_back(transcript_to_junctions.at(t).right);
            }
        } else {
            all_site.right.push_back(transcript_to_junctions.at(t).right);
        }
    }

    return all_site;
}

/*
    returns the value in list that is closest to num
*/
int 
take_closest(const std::vector<int> &list, int num) {
    int output = list[0];

    for (const auto & i : list) {
        if (abs(i - num) < abs(output - num)) {
            output = i;
        }
    }
    return output;
}


std::set<int> 
get_splice_site(
    const std::unordered_map<std::string, Junctions> &transcript_to_junctions,
    const std::vector<std::string> &tr_list)
{
    std::set<int> all_site;

    // add all of the junctions of 
    // everything from tr_list to all_site
    for (const std::string &t : tr_list) {
        for (const auto & junctions : transcript_to_junctions.at(t).junctions) {
            all_site.insert(junctions);
        }
    }

    return all_site;
}


std::vector<int>
find_best_splice_chain(const std::vector<int> &raw_iso, const std::vector<std::vector<int>> &junction_list, int max_dist)
{
    // index, length-1, starting position
    int best_match[3] = {-1, 0, 0};
    for (int i = 0; i < (int)junction_list.size(); i++) {
        std::vector<int> junction = junction_list[i];

        // populate vector i_st with the indices of junction_list[i] entries
        // that match certain criteria
        std::vector<int> iter_start;
        for (int j = 0; j < (int)junction_list[i].size(); j++) {
            if (abs(junction[j] - raw_iso[1]) < max_dist) {
                iter_start.push_back(j);
            }
        }

        if (iter_start.size() == 1) {
            int iter_start_int = iter_start[0];
            int iter = iter_start_int + 1;
            int iter_end = iter_start_int;
            while (iter < (int)junction.size() && iter - iter_start_int + 1 < (int)raw_iso.size()) {
                if (abs(junction[iter] - raw_iso[iter - iter_start_int + 1]) < max_dist) {
                    iter_end = iter++;
                } else {
                    break;
                }
            }

            if (iter_end - iter_start_int >= best_match[1]) {
                best_match[0] = i;
                best_match[1] = iter_end - iter_start_int;
                best_match[2] = iter_start_int;
            }
        }
    }

    if (best_match[0] >= 0 && best_match[1] >= 3) {
        std::vector<int> updated_iso = raw_iso;

        for (int i = best_match[2]; i < best_match[2] + best_match[1] + 1; i++) {
            updated_iso[i - best_match[2] + 1] = junction_list[best_match[0]][i];
        }

        return updated_iso;
    } else {
        return raw_iso;
    }
}


/*
    takes a vector,
    splits it up into a vector of StartEndPairs
    {1, 2, 3, 4, 5} -> {{1, 2}, {3, 4}}
*/
std::vector<exon>
pairwiseStartEndPair (const std::vector<int> &input)
{
    std::vector<exon> output;

    for (int i = 1; i < (int)input.size(); i += 2) {
        StartEndPair new_pair = {input[i-1], input[i]};

        output.push_back(new_pair);
    }

    return output;
}

/* takes two exons, returns the total overlap between them */
int
exon_overlap(const std::vector<exon> &exons1, const std::vector<exon> &exons2)
{
    int total = 0;
    for (const auto & e1 : exons1) {
        for (const auto & e2 : exons2) {
            total += 
                std::max(0, 
                    std::min(e1.end, e2.end) 
                    - std::max(e2.start, e1.start));;
        } 
    }
    return total;
}
int
exon_overlap (const std::vector<int> &exons1, const std::vector<exon> &exons2)
{
    return exon_overlap(pairwiseStartEndPair(exons1), exons2);
}

/*
    checks if s2 is in s1
    only searches for exact matches
*/
bool
check_exon_subset(const std::vector<int> &s1, const std::vector<int> &s2, int MAX_TOLERANCE)
{
    if ((int)s2.size() == 2) { // ignore single exon transcripts 
        return false;
    }

    auto fs_elem = std::find(s1.begin(), s1.end(), s2[1]);

    if (fs_elem == s1.end()) { // ignore if s2[1] is not in s1 
        return false;
    }

    // get the index of the element
    int fs = std::distance(s1.begin(), fs_elem);

    if ((fs == 0) || ((s2[0] - s1[fs - 1]) < -MAX_TOLERANCE)) { // ignore if left is not within s1 
        return false;
    }

    for (int i = 2; i < (int)s2.size() - 1; i++) {
        if (fs + i - 1 > (int)s1.size() - 1) {
            return false;
        }

        if (s1[fs + i - 1] != s2[i]) {
            return false;
        }
    }

    if (fs + (int)s2.size() - 2 > (int)s1.size() - 1) {
        return false;
    }
    if ((s2.back() - s1[fs + (int)s2.size() - 2]) > MAX_TOLERANCE) {
        return false;
    }

    // if all of this is true, return true
    return true;
}

/*
    takes two exon transcripts,
    returns the percentage of coverage between them
*/
float
get_exon_sim_pct(const std::vector<int> &exons1, const std::vector<int> &exons2)
{
    /* add up the pairs of an exon */
    auto
    sum_of_exon = [] (const std::vector<int> &ex)
    {
        return ranges::sumMap<exon, int>(
            pairwiseStartEndPair(ex), 
            [](exon sep) { return sep.end - sep.start; });
    };

    auto e1_len = sum_of_exon(exons1);
    auto e2_len = sum_of_exon(exons2);

    float total = 0;
    for (const auto & pair1 : pairwiseStartEndPair(exons1)) {
        for (const auto & pair2 : pairwiseStartEndPair(exons2)) {
            if ((pair1.end > pair2.start) &&
                (pair1.start < pair2.end)) {
                total += std::min(pair1.end, pair2.end) - std::max(pair1.start, pair2.start);
            } 
        }
    }

    return total / std::max(e1_len, e2_len);
}