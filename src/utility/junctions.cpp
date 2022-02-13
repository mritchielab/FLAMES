/*
  Anything related to junctions are goes in here
*/
#include "junctions.h"

#include <stdio.h>
#include <iostream>
#include <map>
#include <vector>
#include <any>
#include <set>
#include <list>
#include <algorithm>

#include "../classes/GeneBlocks.h"
#include "../classes/StartEndPair.h"


int 
take_closest(std::vector<int> list, int num) {
    /*
        returns the value in list that is closest to num
    */

    int output = list.back();
    list.pop_back();

    for (const auto & i : list) {
        if (abs(i - num) < abs(output - num)) {
            output = i;
        }
    }
    return output;
}


Junctions 
blocks_to_junctions (std::vector<StartEndPair> blocks)
{
    /*
        takes the blocks,
        converts them into a junctions object
    */

    Junctions output;

    // output.left = {blocks.front().start};
    // output.right = {blocks.back().end};
	output.left = blocks.front().start;
	output.right = blocks.back().end;

    if (blocks.size() > 1) {
        for (int i = 1; i < blocks.size(); i++) {
            output.junctions.push_back(blocks[i - 1].end);
            output.junctions.push_back(blocks[i].start);
        }
    }

    return output; 
}


DoubleJunctions 
get_TSS_TES_site
(
    const std::unordered_map<std::string, Junctions> &transcript_to_junctions,
    const std::vector<std::string> &tr_list
)
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


std::set<int> 
get_splice_site 
(
    std::unordered_map<std::string, Junctions> transcript_to_junctions,
    std::vector<std::string> tr_list
)
{
    std::set<int>
    all_site;

    // add all of the junctions of 
    // everything from tr_list to all_site
    for (std::string t : tr_list) {
        for (const auto & junctions : transcript_to_junctions[t].junctions) {
            all_site.insert(junctions);
        }
    }

    return all_site;
}

std::unordered_map<std::string, std::vector<StartEndPair>>
get_gene_flat
(
    std::unordered_map<std::string, std::vector<std::string>>   * gene_to_transcript,
    std::unordered_map<std::string, std::vector<StartEndPair>>  * transcript_to_exon
)
{
    std::unordered_map<std::string, std::vector<StartEndPair>>
    gene_dict = {};

    for (const auto & [gene, transcript] : (*gene_to_transcript)) {
        gene_dict[gene] = (*transcript_to_exon)[transcript[0]];

        if (transcript.size() > 1) {
            // for each entry in transcript
            for (auto t = transcript.begin() + 1; t != transcript.end(); ++t) {
                // and then for each entry in the corresponding exon
                for (auto exon = (*transcript_to_exon)[*t].begin(); exon != (*transcript_to_exon)[*t].end(); ++exon) {
                    bool
                    exon_found = false;
                    // and then for each entry in the gene dictionary
                    for (auto pair = gene_dict[gene].begin(); pair != gene_dict[gene].end(); ++pair) {
                        // disregard ones that have no overlap whatsoever
                        if ((exon->start > pair->end) or (exon->end < pair->start)) {
                            continue;
                        }

                        // otherwise, we've found an exon
                        exon_found = true;
                        // so update the dictionary values accordingly
                        if (exon->end > pair->end) {
                            pair->end = exon->end;
                        }
                        if (exon->start < pair->start) {
                            pair->start = exon->start;
                        }
                    }
                    
                    // if there was no overlapping, just add it as a new exon
                    if (!exon_found) {
                        gene_dict[gene].push_back((*exon));
                    }
                }
            }
        }
        //sort by starting position
        std::sort(gene_dict[gene].begin(), gene_dict[gene].end(), 
            [] (const auto & p1, const auto & p2) {
                return (p1.start < p2.start);
            }
        );
    }

    return gene_dict;
}

void
remove_similar_tr
(
    std::unordered_map<std::string, std::vector<std::string>>   &gene_to_transcript,
    const std::unordered_map<std::string, std::vector<StartEndPair>>  &transcript_to_exon,
    int threshold
)
{
    int
    dup_count = 0, t1, t2;
    
    for (const auto &[gene, transcript] : gene_to_transcript) {
        // ignore anything short
        if (transcript.size() < 2) {
            continue;
        }

        // keep a set of all the duplicates
        std::set<int>
        dup_set = {};

        // iterate through every pair of transcripts t1 and t2 in the gene
        for (t1 = 0; t1 < transcript.size() - 1; ++t1) {
            for (t2 = t1 + 1; t2 < transcript.size(); ++t2) {
            // std::cout << "t2 is " << t2 << " / " << transcript.size() << "\n";
                if (is_exon_similar(transcript_to_exon.at(transcript[t1]), transcript_to_exon.at(transcript[t2]), threshold)) {
                    dup_set.insert(t2);
                    dup_count++;
                }
            }
        }
        
        // then remove the duplicates
        if (dup_count > 0) {
            for (const auto & i : dup_set) {
                gene_to_transcript[gene].erase(transcript.begin() + i);
            }
        }
    }
}

bool
is_exon_similar
(
    const std::vector<StartEndPair> &exon1, 
    const std::vector<StartEndPair> &exon2, 
    int threshold
)
{
    /* 
        takes two exons,
        checks that they are similar enough within a given threshold
    */

    // make sure they're the same size
    if (exon1.size() != exon2.size()) {
        return false;
    }

    // then add up the difference
    int diff = 0;
    for (int i = 0; i < exon1.size(); ++i) {
        diff += std::abs(exon1.at(i).start - exon2.at(i).start) + std::abs(exon1.at(i).end - exon2.at(i).end);
        if (diff > threshold) {
            return false;
        }
    }
    return true;
}

std::unordered_map<std::string, std::vector<GeneBlocks>>
get_gene_blocks
(
    std::unordered_map<std::string, std::vector<StartEndPair>>  * gene_dict,
    std::unordered_map<std::string, std::vector<std::string>>   * chr_to_gene,
    std::unordered_map<std::string, std::vector<std::string>>   * gene_to_transcript
)
{
    std::unordered_map<std::string, std::vector<GeneBlocks>>
    chr_to_blocks = {};

    for (const auto & [chr, genes] : (*chr_to_gene)) {
        chr_to_blocks[chr] = {};

        struct GeneListEntry {
            int start, end;
            std::vector<std::string>
            transcript_list;
            std::string
            gene;
        };
        // extract all the blocks from gene_dict and gene_to_transcript
        std::vector<GeneListEntry>
        gene_list;
        for (const auto & gene : genes) {
            gene_list.push_back({
                (*gene_dict)[gene].front().start, 
                (*gene_dict)[gene].back().end, 
                (*gene_to_transcript)[gene], 
                gene
            });
        }
        // for (const auto & entry : gene_list)
        // sort them by starting position
        std::sort(gene_list.begin(), gene_list.end(), 
            [] (const auto & p1, const auto & p2) {
                return (p1.start < p2.start);
            }
        );

        // add the first item
        chr_to_blocks[chr].push_back(GeneBlocks(
            gene_list[0].start,
            gene_list[0].end,
            gene_list[0].transcript_list,
            gene_list[0].gene
        ));

        if (gene_list.size() > 1) {
            for (auto entry = gene_list.begin(); entry != gene_list.end(); ++entry) {
                if (chr_to_blocks[chr].back().end > entry->start) {
                    // just edit an existing entry
                    chr_to_blocks[chr].back().add_gene(
                        entry->start,
                        entry->end,
                        entry->transcript_list,
                        entry->gene
                    );
                } else {
                    // add a new entry
                    chr_to_blocks[chr].push_back(GeneBlocks(
                        entry->start,
                        entry->end,
                        entry->transcript_list,
                        entry->gene
                    ));
                }
            }
        }
    }
    return chr_to_blocks;
}