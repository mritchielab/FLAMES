#include "Isoforms.h"

#include <map>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <sstream>
#include <iostream>
#include <tuple>
#include <fstream>
#include <iterator>

#include "../utility/junctions.h"
#include "../utility/misc.h"
#include "../utility/utility.h"
#include "GeneBlocks.h"
#include "Config.h"
#include "StartEndPair.h"
#include "Pos.h"

/* 
  add one new isoform to this Isoforms object,
  by either adding it or updating an existing entry that is close to it
*/
void Isoforms::add_isoform(Junctions junctions, bool is_reversed) {
	if (junctions.junctions.size() == 0) { // single exon reads 
		if (this->single_block_dict.size() == 0) {
			this->add_one(junctions, is_reversed);
		} else {
			// there is already stuff in single_block_dict
			// first, check if this key is already in there
			StartEndPair target_key = {junctions.left, junctions.right};

			if (this->single_block_dict.count(target_key) > 0) { // the key is already present
				this->update_one(junctions, target_key, is_reversed);
			} else { // the key is not present, so check for similar keys
				int max_ts_dist = this->parameters.MAX_TS_DIST;
				auto found_similar = std::find_if(single_block_dict.begin(), single_block_dict.end(), 
					[&target_key, &max_ts_dist](std::pair<StartEndPair, std::vector<StartEndPair>> curr) { 
						return ((std::abs(curr.first.start - target_key.start) < max_ts_dist) &&
							(std::abs(curr.first.end - target_key.end) < max_ts_dist)); });

				if (found_similar != single_block_dict.end()) {
					this->update_one(junctions, found_similar->first, is_reversed);
				} else {
					this->add_one(junctions, is_reversed);
				}

				// bool found_similar = false;
				// for (auto const & [current_key, _] : this->single_block_dict) {
				// 	if ((std::abs(current_key.start - target_key.start) < this->parameters.MAX_TS_DIST) &&
				// 			(std::abs(current_key.end - target_key.end) < this->parameters.MAX_TS_DIST)) {
				// 		this->update_one(junctions, current_key, is_reversed);
				// 		found_similar = true;
				// 		break;
				// 	}
				// }

				// if (!found_similar) {
				// 	this->add_one(junctions, is_reversed);
				// }
			}
		}
	} else { // multi-exon reads 
		if (this->junction_dict.size() == 0) {
			this->add_one(junctions, is_reversed);
		} else {
			// there's stuff in junction_dict
			// check to see if there's an exact match
			if (this->junction_dict.count(junctions.junctions) > 0) {
				// the key is not present in the dict yet
				// just add it
				this->update_one(junctions, junctions.junctions, is_reversed);
			} else {
				// the key is not directly present in the dict
				// check for similar keys
				bool found = false;
				for (auto & [current_key, _] : this->junction_dict) {
					if (junctions.junctions.size() == current_key.size()) {
						// they are the same size. see if they are similar
						// instead of comparing
						bool similar = true;
						for (int i = 0; i < current_key.size(); i++) {
							if (std::abs(current_key[i] - junctions.junctions[i]) >= this->parameters.MAX_DIST) {
								// they are too far apart
								similar = false;
								break;
							}
						}

						if (similar) {
							// if they are still similar,
							// then we have found an entry that is close enough
							this->update_one(junctions, current_key, is_reversed);
							found = true;
							break;
						}
					}
				}
				
				if (!found) {
					// we went through and found nothing similar enough
					// so just add it to the dict
					this->add_one(junctions, is_reversed);
				}
			}
		}
	}
};

/* takes a junctions object,
 * adds it to junction_list, junction_dict,
 * left, right, and lr_pairs
 */
void Isoforms::add_one(Junctions junctions, bool strand) {
	int multiplier = strand ? -1 : 1;

  	if (junctions.junctions.size() == 0) { // single-exon
   		StartEndPair key = {junctions.left, junctions.right};

		// this->single_block_dict[key] = this->single_blocks.size();
		// this->single_blocks.push_back({key});
		this->single_block_dict[key] = {key};

		this->strand_counts[std::vector<int> {key.start, key.end}] = {multiplier * this->parameters.STRAND_SPECIFIC};
	} else { // multi-exon reads
		// this->junction_dict[junctions.junctions] = this->junction_list.size();
		// this->junction_list.push_back({junctions.junctions});
		this->junction_dict[junctions.junctions] = std::vector<std::vector<int>>{junctions.junctions};
		this->left.push_back(junctions.left);
		this->right.push_back(junctions.right);
		this->lr_pair[junctions.junctions] = { StartEndPair {junctions.left,junctions.right} };

		this->strand_counts[junctions.junctions] = { multiplier * this->parameters.STRAND_SPECIFIC };
	}
}

/* takes a junctions object and a key to update,
 * updates the key entry in single_blocks, junction_list, left, right, lr_pair
 * and new_strand_counts with the new key
 * This overload is for single-exon reads, where junctions.junctions.size() is assumed to be 0
 */
void Isoforms::update_one(Junctions junctions, StartEndPair key, bool strand) {
    // this->single_blocks[this->single_block_dict[key]].push_back(key);
	single_block_dict.at(key).push_back({junctions.left, junctions.right});

	int multiplier = strand ? -1 : 1;
  	this->strand_counts[{key.start, key.end}].push_back(multiplier * this->parameters.STRAND_SPECIFIC);
}
/* takes a junctions object and a key to update,
 * updates the key entry in single_blocks, junction_list, left, right, lr_pair
 * and new_strand_counts with the new key
 * This overload is for multi-exon reads, where junctions.junctions.size() is assumed to not be 0
 */
void Isoforms::update_one(Junctions junctions, std::vector<int> key, bool strand) {
	// this->junction_list[this->junction_dict[key]].push_back(junctions.junctions);
	this->junction_dict[key].push_back(junctions.junctions);
    this->left.push_back(junctions.left);
    this->right.push_back(junctions.right);

    this->lr_pair[key].push_back({junctions.left, junctions.right});

	int multiplier = strand ? -1 : 1;
  	this->strand_counts[key].push_back(multiplier * this->parameters.STRAND_SPECIFIC);
}

/* returns the length of the Isoforms object 
 */
int Isoforms::size() {
  return (this->junction_dict.size() + this->single_block_dict.size());
}


/* go through junction_dict, single_block_dict, strand_counts_dict and lr_pair
 * and update any entries that need updating
 */
void Isoforms::update_all_splice() {
	// std::unordered_map<std::vector<int>, int>
	std::unordered_map<std::vector<int>, std::vector<std::vector<int>>> 
	junction_tmp = this->junction_dict;
	this->junction_dict.clear();

	std::unordered_map<std::vector<int>, std::vector<int>> 
	strand_counts_tmp = this->strand_counts;
	this->strand_counts.clear();

	std::unordered_map<std::vector<int>, std::vector<StartEndPair>>
	lr_pair_tmp = this->lr_pair;
	this->lr_pair.clear();
	// look through junction_dict and update entries
	// with the most common value in each vector
	for (auto const& [key, junction_vec] : junction_tmp) {
		if (junction_vec.size() >= this->parameters.MIN_SUP_CNT) {
			// produce a new key from the most common
			// junction at each position in the vector of junction vectors
			std::vector<int> newKey = mostCommonEachCell(junction_vec, key.size());
			int newStrandValue = mostCommon<int>(strand_counts_tmp.at(key));

			this->junction_dict[newKey] = junction_vec;
			this->strand_counts[newKey] = {newStrandValue};
			this->lr_pair[newKey] = lr_pair_tmp.at(key);
		}
	}

	std::unordered_map<StartEndPair, std::vector<StartEndPair>> 
	single_block_tmp = this->single_block_dict;
	this->single_block_dict.clear();
	// look through single_block_dict and update entries
	// that pass a min count test 
	// with the most common value in each vector
	for (auto const & [key, blocks] : single_block_tmp) {
		if (blocks.size() >= this->parameters.MIN_SUP_CNT) {
			StartEndPair newKey = mostCommonSEP(blocks); // uses StartEndPair specific mostCommon to find the most common individual start and end values
			int commonCount = mostCommon<int>(strand_counts_tmp.at({key.start, key.end}));
			
			this->single_block_dict[newKey] = single_block_tmp.at(key);
			this->strand_counts[{newKey.start, newKey.end}] = {commonCount};
		}
	}
}

std::vector<std::pair<int, int>> 
Isoforms::filter_site(const std::unordered_map<int, int> &list_counts, float fdr_cutoff) {
	std::vector<std::pair<int, int>> mx (list_counts.begin(), list_counts.end());
	ranges::sort(mx, [](std::pair<int, int> a, std::pair<int, int> b) { 
		if (a.second == b.second) {
	        return a.first < b.first;
	    } 
		return a.second > b.second;
	});

	if (mx[0].second == 1 || mx.size() < 5) {
		return mx;
	}

	int trun_num = std::max(2, int(floor(0.05 * mx.size())));

	// the average of the second values in each pair, cutting off at trunc_num
	float rate = 1.0f / (
		std::accumulate(mx.begin(), mx.end() - trun_num, 0.0f, 
			[] (const auto & p1, const auto & p2) {
				return (p1 + p2.second);}
		) / (mx.size() - trun_num)
	);

	std::vector<float> prob;
	for (const auto & i : mx) {
		prob.emplace_back(rate * exp(-rate * i.second));
	}
	
	// create and populate a cumulative_probability variable
	std::vector<float> cumulative_probability { prob[0] };
	for (int i = 1; i < prob.size(); i++) {
		cumulative_probability.push_back(cumulative_probability[i-1] + prob[i]);
	}

	if (cumulative_probability[0] > fdr_cutoff) {
		return mx;
	} else {
		std::vector<std::pair<int, int>> sliced_mx;

		for (int i = 0; i < cumulative_probability.size(); i++) {
			// finish if we've gone past the cutoff
			if (cumulative_probability[i] > fdr_cutoff) {
				break;
			}
			
			sliced_mx.emplace_back(mx[i]);
		}
		return sliced_mx;
	}
}

std::vector<int> 
Isoforms::insert_dist(std::vector<int> fs, std::vector<int> known_site) {
	std::vector<int> tmp = {fs[0]};

	if (fs.size() > 1) {
		for (int it = 1; it < fs.size(); it++) {
			int clo_p = take_closest(tmp, fs[it]);

			if (std::abs(clo_p - fs[it]) > this->parameters.MAX_TS_DIST / 2) {
				tmp.push_back(fs[it]);
			}
		}
	}
	
	if (known_site.size() > 0) {
		for (auto s : known_site) {
			int clo_p = take_closest(tmp, s);
			
			// add it to tmp if it's far enough away
			if (std::abs(clo_p - s) > this->parameters.MAX_TS_DIST) {
				tmp.push_back(s);
			}
		}
	}

	return tmp;
}

std::vector<int> 
Isoforms::countLR(std::ofstream &out_f, const std::unordered_map<int, int> &counts, std::vector<int> junctionsLR, float fdr_cutoff) {
	// generic implementation for both sides
	std::vector<std::pair<int, int>> fs = filter_site(counts, fdr_cutoff);

	if (fs.size() == 0) {
		return {-99999999};
	} 

	std::vector<int> out_cnt;
	
	// print everything to out_f
	for (const auto & it : fs) {
		output_to_bedgraph(out_f, this->ch, it.first, it.first + 1, it.second);
	}

	// make a vector to store the keys from fs_r
	std::vector<int> fs_keys
		= ranges::map<std::pair<int, int>, int>(fs, [](const auto &i){ return i.first; });

	out_cnt = insert_dist(fs_keys, junctionsLR);
	std::sort(out_cnt.begin(), out_cnt.end());

	return out_cnt;
}

void
Isoforms::output_to_bedgraph(std::ofstream &out_f, std::string ch, int start, int end, int count) {
	out_f 
		<< ch << "\t" 
		<< start << "\t" 
		<< end << "\t" 
		<< count << "\n";
}

/*take a known site, write it to an output file and include it in raw_isoforms and strand_count
 */
void Isoforms::filter_TSS_TES(std::ofstream &out_f, DoubleJunctions known_site, float fdr_cutoff) {
	
	// make counts objects for both left and right
	std::unordered_map<int, int> left_counts = countUnique(this->left);
	std::unordered_map<int, int> right_counts = countUnique(this->right);

	if ((left_counts.size() < this->parameters.MIN_SUP_CNT) || 
		(right_counts.size() < this->parameters.MIN_SUP_CNT)) {
		return;
	} 

	std::vector<int> cnt_l = countLR(out_f, left_counts, known_site.left, fdr_cutoff);
	std::vector<int> cnt_r = countLR(out_f, right_counts, known_site.right, fdr_cutoff);

	for (const auto &[junction, block] : this->lr_pair) {
		if (block.size() == 0) {
			continue;
		}

		std::vector<std::pair<StartEndPair, int>> tmp_pair
			= sortNumberOccurances(block);
		std::vector<StartEndPair> pair_after_filtering;
		std::vector<std::pair<StartEndPair, int>> pair_enrich;

		for (const auto &[p, _] : tmp_pair) {
			StartEndPair cl_p = {
				take_closest(cnt_l, p.start),
				take_closest(cnt_r, p.end)
			};
			
			if ((std::abs(cl_p.start - p.start) < (0.5 * this->parameters.MAX_TS_DIST)) and
				(std::abs(cl_p.end - p.end) < 0.5 * this->parameters.MAX_TS_DIST)) {
				if (pair_after_filtering.size() > 0) {
					// now we want to extract the left and right elements
					std::vector<int> pair_after_filtering_left;
					std::vector<int> pair_after_filtering_right;
					for (const auto & it : pair_after_filtering) {
						pair_after_filtering_left.push_back(it.start);
						pair_after_filtering_right.push_back(it.end);
					}

					int distance = (
						std::abs(take_closest(pair_after_filtering_left, cl_p.start) - cl_p.start) +
						std::abs(take_closest(pair_after_filtering_right, cl_p.end) - cl_p.end)
					);

					if (distance > this->parameters.MAX_TS_DIST) {
						if ((cl_p.start < junction.front()) and cl_p.end > junction.back()) {
							pair_after_filtering.push_back(cl_p);
						}
					}
				} else {
					if ((cl_p.start < junction.front()) and (cl_p.end > junction.back())) {
						pair_after_filtering.push_back(cl_p);
					}
				}

				// but we only want to allow up to MAX_SITE_PER_SLICE combinations
				if (pair_after_filtering.size() >= this->parameters.MAX_SITE_PER_SPLICE) {
					break;
				}
			}
		}

		if (pair_after_filtering.size() == 0) {
			// then search for isoform-specific enrichment
			for (const auto & [p, _] : tmp_pair) {
				// we need to sum up all of the elements in tmp_pair that meet our criteria
				int sum = ranges::sumMap<std::pair<StartEndPair, int>, int>(
					tmp_pair, 
					[&, p=p](auto x) { 
						return std::abs(x.first.start - p.start) + std::abs(x.first.end - p.end) < parameters.MAX_TS_DIST ? x.second : 0; 
					}
				);
				pair_enrich.push_back({p, sum});
			}
			std::sort(
				pair_enrich.begin(), pair_enrich.end(),
				[] (const auto & p1, const auto & p2) {
				// sign is flipped because we are reverse sorting
					return (p1.second < p2.second);
				}
			);

			for (const auto & [p, _] : pair_enrich) {
				if (pair_after_filtering.size() > 0) {
					// now we want to extract the left and right elements again
					std::vector<int> pair_after_filtering_left;
					std::vector<int> pair_after_filtering_right;

					for (const auto & it : pair_after_filtering) {
						pair_after_filtering_left.push_back(it.start);
						pair_after_filtering_right.push_back(it.end);
					}

					int distance = (
						std::abs(take_closest(pair_after_filtering_left, p.start) - p.start) +
						std::abs(take_closest(pair_after_filtering_right, p.end) - p.end)
					);

					if ((distance > this->parameters.MAX_TS_DIST) and 
						(p.start < junction.front()) and 
						(p.end > junction.back())) {
						// append the key of p
						pair_after_filtering.push_back(p);
					}
				} else if ((p.start < junction.front()) and (p.end > junction.back())) {
					// append the key of p
					pair_after_filtering.push_back(p);
				}

				// we only want up to MAX_SITE_PER_SLICE combinations
				if (pair_after_filtering.size() >= this->parameters.MAX_SITE_PER_SPLICE) {
					// so just disreagard any after that point
					break;
				}
			}
		}

		float support_count_total = 0;
		for (const auto & p : pair_after_filtering) {
			// now we want to add the values of keys matching certain criteria
			support_count_total += ranges::sumMap<std::pair<StartEndPair, int>, int>(
				tmp_pair, [&, p=p](auto x) { 
					return std::abs(x.first.start - p.start) + std::abs(x.first.end - p.end) < parameters.MAX_TS_DIST ? x.second : 0; 
				}
			);
		}

		if ((support_count_total / this->lr_pair.size()) <= this->parameters.MIN_SUP_PCT) {
			// not enough support counts - this often happens in the case when the TSS/TES is strongly degraded

			if (pair_enrich.size() == 0) {
				for (const auto & [p, _] : tmp_pair) {
					// add it to pair_enrich
					// we need to sum up all of the elements in tmp_pair that meet our criteria
					int sum = ranges::sumMap<std::pair<StartEndPair, int>, int>(tmp_pair, [&, p=p](auto it) { 
							return std::abs(it.first.start - p.start) + std::abs(it.first.end - p.end) < this->parameters.MAX_TS_DIST ? it.second : 0; 
						}
					);
					pair_enrich.push_back({p, sum});
				}
				std::sort(
					pair_enrich.begin(), pair_enrich.end(),
					[] (const auto & p1, const auto & p2) {
						// sign is flipped because we are reverse sorting
						return (p1.second < p2.second);
					}
				);
			}


			std::vector<int> tmp_ex = {(int)(pair_enrich[0].first.start)};
			for (int i : junction) {
				tmp_ex.push_back(i);
			}
			tmp_ex.push_back((int)(pair_enrich[0].first.end));

			this->raw_isoforms[tmp_ex] = block.size();
			this->strand_counts[tmp_ex] = this->strand_counts[junction];
		} else {
			// now, add filtered TSS/TES to raw_isoforms
			for (const auto & p : pair_after_filtering) {
				std::vector<int> tmp_ex { (int)(p.start) };
				for (int i : junction) {
					tmp_ex.push_back(i);
				}
				tmp_ex.push_back((int)(p.end));

				// add up the values to assign to raw_isoforms
				int sum = ranges::sumMap<std::pair<StartEndPair, int>, int>(tmp_pair, [&, p=p](auto it) {
					return std::abs(it.first.start - p.start) + std::abs(it.first.end - p.end) < this->parameters.MAX_TS_DIST ? it.second : 0;
				});
				this->raw_isoforms[tmp_ex] = sum;
				this->strand_counts[tmp_ex] = this->strand_counts[junction];
			}
		}
	}
}

void Isoforms::update_new_isoform(const std::vector<int> &key, long sup_cnt, const std::string &transcript_id, const std::string &gene_id) {
	this->new_isoforms[key] = Iso { sup_cnt, transcript_id, gene_id };
}

/*
  this function does most of the heavy lifting in the Isoform class
  it does of stuff - honestly my eyes glaze over when I try to understand it all at once
*/
void Isoforms::match_known_annotation 
(
  const std::unordered_map<std::string, Junctions> &transcript_to_junctions,
  const std::unordered_map<std::string, Pos> &transcript_dict,
  const std::unordered_map<std::string, std::vector<StartEndPair>> &gene_dict,
  GeneBlocks one_block,
  std::unordered_map<std::string, std::string> fa_dict
)
{
	std::set<int> splice_site = get_splice_site(transcript_to_junctions, one_block.transcript_list);

	// if this is a single-exon read, we're done
	if (splice_site.size() == 0) {
		return;
	}

	std::vector<std::vector<int>> junction_list;
	std::unordered_map<std::vector<int>, std::string> exons_dict;

	// populate junction_list, exons_list and exons)dict
	for (const auto & i : one_block.transcript_list) {
		const auto &junction = transcript_to_junctions.at(i);
		junction_list.push_back(junction.junctions);

		// add the new entry to exons_list
		// exons list is every entry from junction_list with left and right on either side
		std::vector<int> cur_exon = {junction.left};
		cur_exon.insert(cur_exon.end(), junction.junctions.begin(), junction.junctions.end());
		cur_exon.push_back(junction.right);
		exons_dict[cur_exon] = i;
	}

	DoubleJunctions
	TSS_TES_site = get_TSS_TES_site(transcript_to_junctions, one_block.transcript_list);


	for (const auto & [exon_key, exon_val] : this->single_block_dict) {
		for (const auto & i : one_block.transcript_list) {
			const auto &junction = transcript_to_junctions.at(i);

			int tmp_std = !this->parameters.STRAND_SPECIFIC ? 
				this->strand_counts[{exon_key.start, exon_key.end}][0] : 
				transcript_dict.at(i).strand == '+' ? 1 : -1;

			if ((junction.junctions.size() == 0) && (tmp_std == this->strand_counts[{exon_key.start, exon_key.end}][0])) {
				if ((std::abs(exon_key.start - junction.left) < this->parameters.MAX_TS_DIST) &&
					(std::abs(exon_key.end - junction.right) < this->parameters.MAX_TS_DIST)) {
					
					std::vector<int> known_exons {junction.left, junction.right};
					long new_support_count;
					if (this->known_isoforms.count(known_exons)) { // check if it is already known
						// it's already known, just update the entry
						new_support_count = this->known_isoforms[known_exons].support_count + this->single_block_dict[exon_key].size();
					} else {
						// it's totally new, create a fresh entry for it
						new_support_count = this->single_block_dict[exon_key].size();

						this->strand_counts[known_exons] = !this->parameters.STRAND_SPECIFIC ? 
							std::vector<int> (transcript_dict.at(i).strand == '+' ? 1 : -1) : 
							this->strand_counts[{exon_key.start, exon_key.end}];
					}

					this->known_isoforms[known_exons] = Iso {
						new_support_count,
						i,
						transcript_dict.at(i).parent_id
					};
				}
			}
		}
	}

	// at this point everything is correct and working
	for (auto const & [raw_iso_key, raw_iso_val] : this->raw_isoforms) {
		bool found = false;

		// we need to check if there are any repeated values in raw_iso_key
		// and ignore raw_iso data that has repeated splice sites
		if (ranges::hasDuplicates(raw_iso_key)) {
			continue;
		}

		for (const auto & i : one_block.transcript_list) {
			const auto junction = transcript_to_junctions.at(i);

			int tmp_std = !this->parameters.STRAND_SPECIFIC ? 
				this->strand_counts[raw_iso_key][0] : 
				transcript_dict.at(i).strand == '+' ? 1 : -1;

			// same number of exons, same strand
			if ((raw_iso_key.size() - 2 == junction.junctions.size()) &&
				(tmp_std == this->strand_counts[raw_iso_key][0])) {
				
				bool iso_is_within_max_dist = true;
				for (int j = 0; j < std::min(raw_iso_key.size() - 1, junction.junctions.size()); j++) {
					if (std::abs(raw_iso_key[j+1] - junction.junctions[j]) > this->parameters.MAX_TS_DIST) {
						// we don't want any to be greater than MAX_DIST
						iso_is_within_max_dist = false;
						break;
					}
				}

				if  (iso_is_within_max_dist &&
					(std::abs(raw_iso_key.front() - junction.left) < this->parameters.MAX_TS_DIST) &&
					(std::abs(raw_iso_key.back() - junction.right) < this->parameters.MAX_TS_DIST)) {
					// populate known_exons with transcript_to_junctions
					std::vector<int> known_exons = {junction.left};
					known_exons.insert(known_exons.end(), junction.junctions.begin(), junction.junctions.end());
					known_exons.push_back(junction.right);

					if (this->known_isoforms.count(known_exons)) {
						this->known_isoforms[known_exons] = {
							(long)(std::max(this->known_isoforms[known_exons].support_count, long(raw_iso_val))),
							i,
							transcript_dict.at(i).parent_id
						};
					} else {
						this->known_isoforms[known_exons] = Iso {
							raw_iso_val,
							i,
							transcript_dict.at(i).parent_id
						};

						if (!this->parameters.STRAND_SPECIFIC) {
							// if it's not strand specific protocal, use annotation
							this->strand_counts[known_exons] = { transcript_dict.at(i).strand == '+' ? 1 : -1};
						} else {
							this->strand_counts[known_exons] = this->strand_counts[raw_iso_key];
						}
					}

					found = true;
					break;
				}
			}
		}

		if (!found) {
			std::vector<int> new_exons = 
				find_best_splice_chain(raw_iso_key, junction_list, this->parameters.MAX_SPLICE_MATCH_DIST);
			
			if (!isStrictlyIncreasing(new_exons)) {
				new_exons = raw_iso_key;
			}

			for (int i = 0; i < new_exons.size(); ++i) {
				int a_site = new_exons[i];

				if (i == 0) { // the left-most site
					int closest = take_closest(TSS_TES_site.left, a_site);

					if ((std::abs(closest - a_site) < this->parameters.MAX_TS_DIST) && 
						(closest < raw_iso_key[i + 1])) {
						new_exons[i] = closest;
					} else {
						new_exons[i] = a_site;
					}
				} else if (0 < i < raw_iso_key.size() - 1) { // a site from the middle somewhere
					std::vector<int> splice_site_vec (splice_site.begin(), splice_site.end());
					
					int closest = take_closest(splice_site_vec, a_site);

					if ((std::abs(closest - a_site) < this->parameters.MAX_SPLICE_MATCH_DIST) &&
						(closest > new_exons.back()) &&
						(closest < raw_iso_key[i + 1])) {
						new_exons[i] = closest;
					} else {
						new_exons[i] = a_site;
					}
				} else { // then it must be the right-most site
					int closest = take_closest(TSS_TES_site.right, a_site);
					
					if ((std::abs(closest - a_site) < this->parameters.MAX_TS_DIST) &&
						(closest > new_exons.back())) {
						new_exons[i] = closest;
					} else {
						new_exons[i] = a_site;
					}
				}
			}

			if (isStrictlyIncreasing(new_exons)) {
				if (this->new_isoforms.count(new_exons) == 0) {
					this->new_isoforms[new_exons] = Iso {
						this->raw_isoforms[raw_iso_key], 
						"", 
						""
					};
					this->strand_counts[new_exons] = this->strand_counts[raw_iso_key];
				} else {
					this->new_isoforms[new_exons] = Iso {
						this->new_isoforms[new_exons].support_count + raw_iso_val, 
						"", 
						""
					};
				}
			} else {
				// exon chain not in ascending order
			}
		}
	}
	
	// remove incomplete transcript (due to 3' bias)
	if (this->parameters.REMOVE_INCOMP_READS) {
		std::set<std::vector<int>> delete_key;
		auto deleteKeysFromNewIsoforms = [this, &delete_key]() { 
			for (const auto & key : delete_key) {
				this->new_isoforms.erase(key);
			}
		};

		for (const auto & [new_isoform_key, new_isoform_val] : this->new_isoforms) {
			// 1. match to known isoform detected and remove new isoform if within known isoform
			if (this->known_isoforms.count(new_isoform_key)) {
				delete_key.insert(new_isoform_key);
				continue;
			}

			for (const auto & [known_isoform_key, known_isoform_val] : this->known_isoforms) {
				int position_modifier = -1;
				if ((new_isoform_key.size() < known_isoform_key.size()) &&
					(if_exon_contains(known_isoform_key, std::vector<int>(new_isoform_key.begin()+2, new_isoform_key.end()), 1))) {
					position_modifier = new_isoform_key[0];
				} else if ((new_isoform_key.size() < known_isoform_key.size()) &&
						(if_exon_contains(known_isoform_key, std::vector<int>(new_isoform_key.begin(), new_isoform_key.end() - 2), 1))) {
					position_modifier = new_isoform_key.back();
				}

				if (position_modifier != -1) {
					std::vector<char> s_l = std::vector<char>( // get a slice of 15 ints just before new_isoform_key.back()
						fa_dict[this->ch].begin() + (position_modifier - 15), 
						fa_dict[this->ch].begin() + (position_modifier)
					);

					std::vector<char> s_r = std::vector<char>( // get a slice of 15 ints just after new_isoform_key.back()
						fa_dict[this->ch].begin() + position_modifier,
						fa_dict[this->ch].begin() + position_modifier + 15
					);
					
					if ((std::count(s_l.begin(), s_l.end(), 'T') > 10) ||
						(std::count(s_l.begin(), s_l.end(), 'A') > 10) ||
						(std::count(s_r.begin(), s_r.end(), 'T') > 10) ||
						(std::count(s_r.begin(), s_r.end(), 'A') > 10)) {
						// remove keys with too many Ts or As at the beginning or end
						delete_key.insert(new_isoform_key);
						break;
					}
				}

				// test: above is save

				if ((new_isoform_key.size() < known_isoform_key.size()) && if_exon_contains(known_isoform_key, new_isoform_key, this->parameters.MAX_TS_DIST)) {
					float sim_pct_sq = pow(get_exon_sim_pct(new_isoform_key, known_isoform_key), 2);
					if (new_isoform_val.support_count < (1 + sim_pct_sq * this->parameters.REMOVE_INCOMP_READS) * known_isoform_val.support_count) {
						// remove keys with too much similarity in coverage
						delete_key.insert(new_isoform_key);
						break;
					}
				}

				// add to delete keys if they are the same
				if (new_isoform_key == known_isoform_key) {
					delete_key.insert(new_isoform_key);
					break;
				}
			}
		}

		// entirely remove the keys in delete_key_set from new_isoforms 
		deleteKeysFromNewIsoforms();
		// reset delete key
		delete_key.clear();

		// 2. match to known isoform in annotation, and remove new_isoform if it's very similar to annotation
		for (const auto & [new_isoform_key, new_isoform_val] : this->new_isoforms) {
			for (const auto & [known_isoform_key, known_isoform_va] : exons_dict) {
				if (!this->known_isoforms.count(known_isoform_key)) {
					if (new_isoform_key.size() < known_isoform_key.size() &&
						if_exon_contains(known_isoform_key, new_isoform_key, this->parameters.MAX_TS_DIST)) {
						auto sim_pct = get_exon_sim_pct(new_isoform_key, known_isoform_key);

						if (sim_pct > 0.95) {
							delete_key.insert(new_isoform_key);
							this->known_isoforms[known_isoform_key] = Iso {
								this->new_isoforms[new_isoform_key].support_count,
								exons_dict[known_isoform_key],
								transcript_dict.at(exons_dict[known_isoform_key]).parent_id
							};

							this->strand_counts[known_isoform_key] = {
								transcript_dict.at(exons_dict[known_isoform_key]).strand == '+' ? 1 :  -1};
						}
					} else if (new_isoform_key == known_isoform_key) {
						delete_key.insert(new_isoform_key);
						this->known_isoforms[known_isoform_key] = Iso {
							this->new_isoforms[new_isoform_key].support_count,
							exons_dict[known_isoform_key],
							transcript_dict.at(exons_dict[known_isoform_key]).parent_id
						};

						this->strand_counts[known_isoform_key] = {
							transcript_dict.at(exons_dict[known_isoform_key]).strand == '+' ? 1 : -1};
					}
				}
			}
		}
		
		// entirely remove the keys in delete_key_set from new_isoforms 
		deleteKeysFromNewIsoforms();

		// 3. Match among new isoform
		if (this->new_isoforms.size() > 1) {
			// reset delete key
			delete_key.clear();

			for (auto it = this->new_isoforms.begin(); std::next(it, 1) != new_isoforms.end(); ++it) {
				for (auto jt = std::next(it, 1); jt != new_isoforms.end(); ++jt) {
					// remove isoform from new isoforms if similar enough to others in the record
					auto removeSimPctSqReads = [&](const auto &key1, const auto &key2) {
						if (key1->first.size() < key2->first.size() && if_exon_contains(key2->first, key1->first, this->parameters.MAX_TS_DIST)) {
							auto sim_pct_sq = pow(get_exon_sim_pct(key2->first, key1->first), 2);
							
							if (this->new_isoforms[key1->first].support_count < sim_pct_sq * this->parameters.REMOVE_INCOMP_READS * key1->second.support_count) {
								delete_key.insert(key1->first);
							}
						}
					};

					removeSimPctSqReads(it, jt);
					removeSimPctSqReads(jt, it);
				}
			}

			// entirely remove the keys in delete_key from new_isoforms
			deleteKeysFromNewIsoforms();
		}
	}

	// match to gene
	std::set<std::vector<int>> delete_key;
	auto deleteKeysFromNewIsoforms = [this, &delete_key]() { 
		for (const auto & key : delete_key) {
			this->new_isoforms.erase(key);
		}
	};
	std::map<std::vector<int>, Iso> update_iso_dict;

	if (this->new_isoforms.size() > 0) {
		for (auto & [isoform_key, isoform_val] : this->new_isoforms) {
			// calculate a sum of the pairs of the isoform
			int iso_len = ranges::sumMap<StartEndPair, int>(pairwise(isoform_key), [](const StartEndPair &i){ 
				return i.end - i.start; 
			});
			
			// make tmp to store exon overlap and gene name, and populate it
			std::vector<std::pair<int, std::string>> tmp;
			for (const auto & [ge, tr] : one_block.gene_to_transcript) {
				tmp.push_back( {exon_overlap(isoform_key, gene_dict.at(ge)), ge} );
			}

			if (this->parameters.STRAND_SPECIFIC) { // if it has strand-specific protocol, use read
				char stnd = this->strand_counts[isoform_key][0] == 1 ? '+' : '-';

				// update tmp, only including values which match the sign
				tmp = ranges::filter<std::pair<int, std::string>>(tmp, [&](auto it){ 
					return transcript_dict.at(one_block.gene_to_transcript[it.second][0]).strand == stnd; 
				});

				// if there's nothing left, just skip
				if (tmp.size() == 0) {
					continue;
				}
			}

			// come back later and check that this works
			std::sort(tmp.begin(), tmp.end(),
				[] (const auto & p1, const auto & p2) {
					return (p1.first > p2.first);
				}
			);

			if (tmp[0].first > 0) {
				if (isoform_key.front() >= gene_dict.at(tmp[0].second).front().start &&
					isoform_key.back() <= gene_dict.at(tmp[0].second).back().end) {
					isoform_val = Iso {isoform_val.support_count, "", tmp[0].second};
				} else if ((exon_overlap(std::vector<int>(isoform_key.begin() + 2, isoform_key.begin() + 4), gene_dict.at(tmp[0].second)) > 0) &&
						(exon_overlap(std::vector<int>(isoform_key.end() - 4, isoform_key.end() - 2), gene_dict.at(tmp[0].second)) > 0)) {
					auto iso_key_copy = isoform_key;
					auto iso_val_copy = isoform_val;
					// function template
					auto updateIsoDictIfExonOverlap = [&](int index1, int index2, const std::vector<int> &exonKey1, const std::vector<int> &exonKey2) {
						bool completed = false;
						if ((index1 - index2 <= this->parameters.MIN_FL_EXON_LEN) && (exon_overlap(exonKey1, gene_dict.at(tmp[0].second)) == 0)) {
							std::string ba = ranges::slice(fa_dict[this->ch], index2, index1);

							if (ranges::count(ba, 'T') > 0.7 * ba.size() || ranges::count(ba, 'A') > 0.7 * ba.size()) {
								delete_key.insert(iso_key_copy);
								update_iso_dict[exonKey2] = Iso {iso_val_copy.support_count, "", tmp[0].second};
							}
							completed = true;
						} 
						return completed;
					};

					if (!updateIsoDictIfExonOverlap(isoform_key[1], isoform_key[0], 
							ranges::slice(isoform_key, 0, 2), 
							ranges::slice(isoform_key, 2, isoform_key.size()))) {
						if (!updateIsoDictIfExonOverlap(isoform_key.rbegin()[0], isoform_key.rbegin()[1], 
								ranges::slice(isoform_key, isoform_key.size()-2, isoform_key.size()), 
								ranges::slice(isoform_key, 0, isoform_key.size()-2))) {
							
							isoform_val = Iso {isoform_val.support_count, "", tmp[0].second};
						}
					}
					
				} else {
					if (tmp[0].first > 0.8 * iso_len) {
						if ((isoform_key[1] - isoform_key[0] < this->parameters.MIN_FL_EXON_LEN) || 
							(isoform_key.rbegin()[0] - isoform_key.rbegin()[1] < this->parameters.MIN_FL_EXON_LEN)) {
							continue; // alignment artifact
						}
						// might be real eRNA
						update_new_isoform(isoform_key, isoform_val.support_count, "", tmp[0].second);
					} else {
						// get the total of elements in isoform_key that are also in splice_count
						int total = ranges::sumMap<int, int>(isoform_key, [&splice_site](int it) { 
							return splice_site.count(it) > 0; 
						});

						if (total > 4) {
							// more than 4 splice site match
							update_new_isoform(isoform_key, isoform_val.support_count, "", tmp[0].second);
						}
					}
				}
			}
		}

		// remove all the delete_key entries from new_isoforms
		deleteKeysFromNewIsoforms();

		for (const auto & [iso_key, iso_val] : this->new_isoforms) {
			// skip if there are no matching genes
			if (iso_val.gene_id == "") {
				continue;
			}

			// add iso_val.gene_id to ge_dict if it's not already there
			if (this->ge_dict.count(iso_val.gene_id) == 0) {
				this->ge_dict[iso_val.gene_id] = {};
			}
			
			// update the entry with iso_key (automatically makes a new vector if no matching key)
			this->ge_dict[iso_val.gene_id].push_back(iso_key);

			if (!this->parameters.STRAND_SPECIFIC) {
				if (this->strand_counts[iso_key][0] == 0) {
					this->strand_counts[iso_key] = {
						transcript_dict.at(one_block.gene_to_transcript[iso_val.gene_id][0]).strand == '+' ? 1 :  -1};
				}
			}
		}
	}
	
	for (const auto & [iso_key, iso_val] : this->known_isoforms) {
		// update the value with the new iso_key
		this->ge_dict[iso_val.gene_id].push_back(iso_key);
	}
}

std::string
Isoforms::isoform_to_gff3(float isoform_pct=-1) {
	std::vector<std::string> gff_rec = {};

	std::map<std::string, int> transcript_id_dict;

	if (this->new_isoforms.size() + this->known_isoforms.size() == 0) {
		return "";
	}

	for (const auto & [g_key, g_val] : this->ge_dict) {
		std::vector<std::string> gff_tmp = {};

		// add up all the support_count from all the genes in g_val
		int total_cnt = 0;
		for (const auto & e : g_val) {
			if (this->new_isoforms.count(e) > 0) {
				total_cnt += this->new_isoforms[e].support_count;
			}
			if (this->known_isoforms.count(e) > 0) {
				total_cnt += this->known_isoforms[e].support_count;
			}
		}

		std::vector<int> firsts, lasts;
		for (const auto & i : g_val) {
			if ((this->known_isoforms.count(i) > 0) || 
				(this->new_isoforms[i].support_count > isoform_pct * total_cnt)) {
				firsts.push_back(i.front());
				lasts.push_back(i.back());
			}
		}
		
		std::stringstream gff_entry;
		gff_entry 
			<< this->ch << '\t' // _ch
			<< "." << '\t' // _sr
			<< "gene" << '\t' // _ty
			<< *(std::min_element(firsts.begin(), firsts.end())) + 1 << '\t' // _st
			<< *(std::max_element(lasts.begin(), lasts.end())) << '\t' // _en
			<< "." << '\t' // _sc
			<< (this->strand_counts[g_val[0]][0] == 1 ? '+' : '-') << '\t' // _stnd
			<< "." << "\t" // _ph
			<< "ID=gene:" << g_key << ";gene_id=" << g_key << ";support_count=" << total_cnt; // _attr
		
		gff_tmp.push_back(gff_entry.str());

		for (const auto & exons : g_val) {
			if (this->new_isoforms.count(exons) && this->known_isoforms.count(exons)) {
				Rcpp::Rcout << "BOTH in new and known\n";
			}

			std::string source;
			long support_count;
			std::stringstream tp_id;
			if (this->new_isoforms.count(exons)) {
				source = "new";
				support_count = this->new_isoforms[exons].support_count;

				if (0.0f < isoform_pct && isoform_pct < 1.0f && support_count < (isoform_pct * total_cnt)) {
					continue;
				}
				// store the key, the min and the max
				tp_id << g_key << "_" << exons[0]+1 << "_" << exons.back();

				if (transcript_id_dict.count(tp_id.str()) == 0) {
					// if it's not yet there, add tp_id to the dict
					transcript_id_dict[tp_id.str()] = 1;
					// update tp_id
					tp_id << "_1";
				} else {
					// otherwise, just update it
					transcript_id_dict[tp_id.str()] += 1;
					// and then update tp_id
					auto tmp = transcript_id_dict[tp_id.str()];
					tp_id << "_" << tmp;
				}
			} else {
				source = "known";
				support_count = this->known_isoforms[exons].support_count;
				tp_id << this->known_isoforms[exons].transcript_id;
			}

			if (support_count < this->parameters.MIN_SUP_CNT) {
				continue;
			}

			std::stringstream gff_entry;
			gff_entry 
				<< this->ch << '\t' // _ch // chromosome
				<< source << '\t' // _sr // source
				<< "transcript" << '\t' // _ty // 
				<< exons[0] + 1 << '\t' // _st // start
				<< exons.back() << '\t' // _en // end
				<< "." << '\t' // _sc 
				<< (this->strand_counts[exons][0] == 1 ? '+' : '-') << '\t' // _stnd // strand
				<< "." << "\t" // _ph
				<< "ID=transcript:" << tp_id.str() << ";transcript_id=" << tp_id.str() << ";Parent=gene:" << g_key << ";support_count=" << support_count << ";source=" << source; // _attr
			gff_tmp.push_back(gff_entry.str());

			int exon_idx = 1;
			for (int i = 0; i < exons.size(); i+=2) {
				std::stringstream gff_entry;
				gff_entry 
					<< this->ch << '\t' // _ch // chromosome
					<< source << '\t' // _sr // source
					<< "exon" << '\t' // _ty // 
					<< exons[i] + 1 << '\t' // _st // start
					<< exons[i+1] << '\t' // _en // end
					<< "." << '\t' // _sc 
					<< (this->strand_counts[exons][0] == 1 ? '+' : '-') << '\t' // _stnd // strand
					<< "." << "\t" // _ph
					<< "exon_id=exon:" << exons[i]+1 << "_" << exons[i+1] << ";Parent=transcript:" << tp_id.str() << ";rank=" << exon_idx; // _attr
				gff_tmp.push_back(gff_entry.str());
				exon_idx++;
			}
		}

		if (gff_tmp.size() > 2) {
			for (const auto & gff : gff_tmp) {
				gff_rec.push_back(gff);
			}
		}
	}

	if (gff_rec.size() > 0) {
		std::stringstream output;
		for (const auto & gff : gff_rec) {
			output << gff << "\n";
		}
		return output.str();
	}

	return "";
}
