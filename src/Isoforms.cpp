#include "Isoforms.h"

/* 
  add one new isoform to this Isoforms object,
  by either adding it or updating an existing entry that is close to it
*/
void Isoforms::add_isoform(Junctions junctions, bool is_reversed)
{
  if (junctions.junctions.size() == 0) { // single exon reads 
    if (this->single_block_dict.size() == 0) {
      // this will be the first thing in single_block_dict
      // so, nothing to check
      this->add_one(junctions, is_reversed);
    } else {
      // there is already stuff in single_block_dict
      // first, check if this key is already in there
      StartEndPair
      target_key = {junctions.left[0], junctions.right[0]};

      if (this->single_block_dict.count(target_key) > 0) {
        // this time the key must be converted to vector
        // we really need to fix this
        std::vector<int> 
        target_key_vector = {target_key.start, target_key.end};

        // the key is already present - so just update it
        this->update_one(junctions, target_key_vector, is_reversed);
      } else {
        // the key is not present in the dict yet
        // check if there's anything very close to the key
        bool found = false;
        for (auto const & [current_key, current_value] : this->single_block_dict) {
          if ((abs(current_key.start - target_key.start) < this->parameters.MAX_TS_DIST) and
              (abs(current_key.end - target_key.end) < this->parameters.MAX_TS_DIST)) {
            // we've found one that's very similar
            // just update it
            this->update_one(junctions, {current_key.start, current_key.end}, is_reversed);
            found = true;
            break;
          }
        }

        if (!found) {
          // we didn't find any similar ones
          // so just add it
          this->add_one(junctions, is_reversed);
        }
      }
    }
  } else { // multi-exon reads 
    if (this->junction_dict.size() == 0) {
      // there's nothing in the junction dict yet
      // just add it
      this->add_one(junctions, is_reversed);
    } else {
      // there's stuff in junction_dict
      // check to see if there's an exact match

      if (this->junction_dict.count(junctions.junctions) == 0) {
        // the key is not present in the dict yet
        // just add it
        this->add_one(junctions, is_reversed);
      } else {
        // the key is not directly present in the dict
        // check for similar keys

        bool found = false;

        for (auto & [current_key, current_value] : this->junction_dict) {
          if (junctions.junctions.size() == current_key.size()) {
            // they are the same size. see if they are similar
            
            bool similar = true;
            for (int i = 0; i < current_key.size(); i++) {
              if (abs(junctions.junctions[i] - current_key[i]) > this->parameters.MAX_DIST) {
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
void Isoforms::add_one(Junctions junctions, bool strand)
{
  std::map<std::vector<int>, std::vector<int>>
  new_strand_counts;

  if (junctions.junctions.size() == 0) { // single-exon
    StartEndPair key = 
    {junctions.left[0], junctions.right[0]};

    this->single_block_dict[key] = this->single_blocks.size();
    this->single_blocks.push_back({key});
    this->new_strand_counts[std::vector<int> {key.start, key.end}] = {};

    int multiplier = 1;
    if (strand) {
      multiplier = -1;
    }
    this->new_strand_counts[std::vector<int> {key.start, key.end}].push_back(multiplier * this->parameters.STRAND_SPECIFIC);
  } else { // multi-exon reads
    this->junction_dict[junctions.junctions] = this->junction_list.size();
    this->junction_list.push_back({junctions.junctions});
    this->left.push_back(junctions.left[0]);
    this->right.push_back(junctions.right[0]);
    this->lr_pair[junctions.junctions] = {};

    this->new_strand_counts[junctions.junctions] = {};
    int multiplier = 1;
    if (strand) {
      multiplier = -1;
    }
    this->new_strand_counts[junctions.junctions].push_back(multiplier * this->parameters.STRAND_SPECIFIC);
  }
}

/* takes a junctions object and a key to update,
 * updates the key entry in single_blocks, junction_list, left, right, lr_pair
 * and new_strand_counts with the new key
 */
void Isoforms::update_one(Junctions junctions, std::vector<int> key, bool strand)
{
  if (junctions.junctions.size()==0) { // single-exon reads
    std::pair<int, int>
    new_element = {junctions.left[0], junctions.right[0]};

    StartEndPair
    key_pair = {key[0], key[1]};

    this->single_blocks[this->single_block_dict[key_pair]].push_back(key_pair);
  } else { // multi-exon reads 
    this->junction_list[this->junction_dict[key]].push_back(junctions.junctions);
    this->left.push_back(junctions.left[0]);
    this->right.push_back(junctions.right[0]);
    this->lr_pair[key].push_back({junctions.left[0], junctions.right[0]});
  }
  int multiplier = 1;
  if (this->parameters.STRAND_SPECIFIC) {
    multiplier = -1;
  }
  this->new_strand_counts[key].push_back(multiplier * this->parameters.STRAND_SPECIFIC);
}

/* returns the length of the Isoforms object 
 */
int Isoforms::len()
{
  return (this->junction_dict.size() + this->single_block_dict.size());
}


/* go through junction_dict, single_block_dict, strand_counts_dict and lr_pair
 * and update any entries that need updating
 */
void Isoforms::update_all_splice()
{

  std::unordered_map<std::vector<int>, int>
  junction_tmp;
  std::unordered_map<StartEndPair, int>
  single_block_tmp;
  std::unordered_map<std::vector<int>, int>
  strand_counts_tmp;
  std::unordered_map<std::vector<int>, std::vector<StartEndPair>>
  lr_pair_tmp;

  // look through junction_dict
  for (auto const& [key, val] : this->junction_dict) {
    if (this->junction_list[val].size() >= this->parameters.STRAND_SPECIFIC) {
      // variable to store the counts
      std::unordered_map<std::vector<int>, int>
      junction_list_counts;

      // populate the counts container
      for (auto i : this->junction_list[val]) {
        junction_list_counts[i]++;
      }
      
      // pull out the key for the 
      // most common value in the counts map
      auto new_key = (*std::max_element(
        junction_list_counts.cbegin(), junction_list_counts.cend()
      )).first;

      junction_tmp[new_key] = this->junction_dict[key];
      
      // the name is a bit silly, but this is to count the occurrences of keys in strand_counts
      std::unordered_map<int, int>
      strand_counts_counts;
      for (auto i : this->new_strand_counts[key]) {
        strand_counts_counts[i]++;
      }
      // get the key for the most common
      // value in max_strand_counts
      auto max_strand_counts = (*std::max_element(
        strand_counts_counts.cbegin(), strand_counts_counts.cend()
      )).first;

      strand_counts_tmp[new_key] = max_strand_counts;
      lr_pair_tmp[new_key] = this->lr_pair[key];
    }
  }

  // look through single_block_dict
  for (auto const& [key, val] : this->single_block_dict) {
    if (this->single_blocks[val].size() >= this->parameters.STRAND_SPECIFIC) {
      // make and populate a counter for single_blocks
      std::unordered_map<StartEndPair, int>
      single_blocks_counts;
      for (auto i : this->single_blocks[val]) {
        single_blocks_counts[i]++;
      }

      auto new_key = (*std::max_element(
        single_blocks_counts.cbegin(), single_blocks_counts.cend()
      )).first;

      single_block_tmp[new_key] = this->single_block_dict[key];
      std::unordered_map<int, int>
      strand_counts_counts;
      for (auto i : this->new_strand_counts[{key.start, key.end}]) {
        strand_counts_counts[i]++;
      }

      auto max_strand_counts = (*std::max_element(
        strand_counts_counts.cbegin(), strand_counts_counts.cend()
      )).first;

      strand_counts_tmp[{new_key.start, new_key.end}] = max_strand_counts;
    }
  }

  this->single_block_dict = single_block_tmp;
  this->junction_dict = junction_tmp;
  this->strand_counts = strand_counts_tmp;
  this->lr_pair = lr_pair_tmp;
}

/* take a known site, write it to an output file and include it in raw_isoforms and strand_count
 */
void Isoforms::filter_TSS_TES(std::ofstream * out_f, Junctions known_site, float fdr_cutoff)
{
  std::string bedgraph_fmt = "{_ch}\t{_st}\t{_en}\t{_sc}\n";

  auto filter_site = [fdr_cutoff] (std::unordered_map<int, int> list_counts) {
    // to sort the map, we need it as a vector first
    std::vector<std::pair<int, int>>
    mx;
    for (const auto & i : list_counts) {
      mx.emplace_back(i);
    }

    // sort it
    std::sort(mx.begin(), mx.end(),
      [] (const auto & p1, const auto & p2) {
        return (p1.second > p2.second);
      }
    );

    if (mx[0].second == 1 || mx.size() < 5) {
      return mx;
    }

    int trun_num = std::max(2, int(floor(0.05 * mx.size())));

    float rate = 1/(
      // the average of the second values in each pair, cutting off at trunc_num
      std::accumulate(mx.begin(), mx.end() - trun_num, 0, 
        [] (const auto & p1, const auto & p2) {
          return (p1 + p2.second);
        }
      ) / mx.size()
    );

    std::vector<int> 
    prob;

    for (const auto & i : mx) {
      prob.emplace_back(rate * exp(-rate * i.second));
    }
    
    // create and populate a cumulative_probability variable
    std::vector<int>
    cumulative_probability;
    int running_total;
    for (const auto & it : prob) {
      running_total += it;
      cumulative_probability.push_back(running_total);
    }

    if (cumulative_probability[0] > fdr_cutoff) {
      return mx;
    } else {
      std::vector<std::pair<int, int>>
      sliced_mx;

      for (const auto & i : cumulative_probability) {
        // finish if we've gone past the cutoff
        if (i > fdr_cutoff) {
          break;
        }
        
        sliced_mx.emplace_back(mx[i]);
      }
      return sliced_mx;
    }
  };

  auto insert_dist = [this] (std::vector<int> fs, std::vector<int> known_site) 
  {
    std::vector<int>
    tmp = {fs[0]};

    if (fs.size() > 1) {
      for (int it = 1; it != fs.size(); it++) {
        int
        clo_p = take_closest(tmp, fs[it]);

        if (abs(clo_p - fs[it]) > this->parameters.MAX_TS_DIST / 2) {
          tmp.push_back(fs[it]);
        }
      }
    }

    // if known_site is not None
    if (known_site.size() > 0) {
      for (auto s : known_site) {
        int clo_p = take_closest(tmp, s);
        
        // add it to tmp if it's far enough away
        if (abs(clo_p - s) > this->parameters.MAX_TS_DIST) {
          tmp.push_back(s);
        }
      }
    }

    return tmp;
  };
  
  // make counts objects for both left and right
  std::unordered_map<int, int>
  left_counts;
  for (const auto & i : this->left) {
    left_counts[i]++;
  }

  std::vector<int>
  cnt_l;

  std::unordered_map<int, int>
  right_counts;
  for (const auto & i : this->right) {
    right_counts[i]++;
  }

  std::vector<int>
  cnt_r;

  if ((left_counts.size() < this->parameters.STRAND_SPECIFIC) || 
      (right_counts.size() < this->parameters.STRAND_SPECIFIC)) {
    return;
  } else {
    // left
    std::vector<std::pair<int, int>>
    fs_l = filter_site(left_counts);

    if (fs_l.size() == 0) {
      cnt_l = {-99999999};
    } else {
      // print everything to out_f
      for (const auto & it : fs_l) {
        (*out_f) << this->ch << "\t" 
              << it.first << "\t" 
              << it.first + 1 << "\t" 
              << it.second << "\n";
      }

      // make a vector to store the keys from fs_l
      std::vector<int> 
      fs_l_keys;
      for (const auto & i : fs_l) {
        fs_l_keys.push_back(i.first);
      }

      cnt_l = insert_dist(fs_l_keys, known_site.left);
      std::sort(cnt_l.begin(), cnt_l.end());
    }

    // right
    std::vector<std::pair<int, int>>
    fs_r = filter_site(right_counts);

    if (fs_r.size() == 0) {
      cnt_r = {-99999999};
    } else {
      // print everything to out_f
      for (const auto & it : fs_r) {
        (*out_f) << this->ch << "\t"
              << it.first << "\t"
              << it.first + 1 << "\t"
              << it.second << "\n";
      }

      // make a vector to store the keys from fs_r
      std::vector<int>
      fs_r_keys;
      for (const auto & i : fs_r) {
        fs_r_keys.push_back(i.first);
      }

      cnt_r = insert_dist(fs_r_keys, known_site.right);
      std::sort(cnt_r.begin(), cnt_r.end());
    }
  }

  for (const auto & pair : lr_pair) {
    // first we need a counts map, which we populate
    std::unordered_map<StartEndPair, int>
    pair_counts;
    for (const auto & i : pair.second) {
      pair_counts[i]++;
    }

    // then to sort it, we need to convert to a vector
    std::vector<std::pair<StartEndPair, int>>
    tmp_pair;
    for (const auto & i : pair_counts) {
      tmp_pair.push_back(i);
    }

    std::sort(tmp_pair.begin(), tmp_pair.end(),
      [] (const auto & p1, const auto & p2) {
        return (p1.second > p2.second);
      }
    );

    std::vector<StartEndPair>
    pair_after_filtering;

    // lines like these make me think, 
    // maybe i've not thought this architecture out well enough...
    std::vector<std::pair<StartEndPair, int>>
    pair_enrich;

    for (const auto & p : tmp_pair) {
      StartEndPair
      cl_p = {
        take_closest(cnt_l, p.first.start),
        take_closest(cnt_r, p.first.end)
      };

      if ((abs(cl_p.start - p.first.start) < (0.5 * this->parameters.MAX_TS_DIST)) and
          (abs(cl_p.end - p.first.end) < 0.5 * this->parameters.MAX_TS_DIST)) {
        if (pair_after_filtering.size() > 0) {
          // line 650 of sc_longread.py

          // now we want to extract the left and right elements
          std::vector<int>
          pair_after_filtering_left;
          std::vector<int>
          pair_after_filtering_right;

          for (const auto & it : pair_after_filtering) {
            pair_after_filtering_left.push_back(it.start);
            pair_after_filtering_right.push_back(it.end);
          }

          int distance = (
            abs(take_closest(pair_after_filtering_left, cl_p.start) - cl_p.start) +
            abs(take_closest(pair_after_filtering_right, cl_p.end) - cl_p.end)
          );

          if (distance > this->parameters.MAX_TS_DIST) {
            if ((cl_p.start < pair.first.front()) and cl_p.end > pair.first.back()) {
              pair_after_filtering.push_back(cl_p);
            }
          }
        } else {
          if ((cl_p.start < pair.first.front()) and (cl_p.end > pair.first.back())) {
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
      for (const auto & p : tmp_pair) {
        // we need to sum up all of the elements in tmp_pair that meet our criteria
        int sum = 0;
        for (const auto & it : tmp_pair) {
          if ((abs(it.first.start - p.first.start) + abs(it.first.end - p.first.end)) < this->parameters.MAX_TS_DIST) {
            sum += it.second;
          }
        }
        pair_enrich.push_back({p.first, sum});
      }
      std::sort(
        pair_enrich.begin(), pair_enrich.end(),
        [] (const auto & p1, const auto & p2) {
          // sign is flipped because we are reverse sorting
          return (p1.second < p2.second);
        }
      );

      for (const auto & p : pair_enrich) {
        if (pair_after_filtering.size() > 0) {
          // now we want to extract the left and right elements again
          std::vector<int>
          pair_after_filtering_left;
          std::vector<int>
          pair_after_filtering_right;

          for (const auto & it : pair_after_filtering) {
            pair_after_filtering_left.push_back(it.start);
            pair_after_filtering_right.push_back(it.end);
          }

          int distance = (
            abs(take_closest(pair_after_filtering_left, p.first.start) - p.first.start) +
            abs(take_closest(pair_after_filtering_right, p.first.end) - p.first.end)
          );

          if ((distance > this->parameters.MAX_TS_DIST) and 
              (p.first.start < pair.first.front()) and 
              (p.first.end > pair.first.back())) {
            // append the key of p
            pair_after_filtering.push_back(p.first);
          }
        } else if ((p.first.start < pair.first.front()) and (p.first.end > pair.first.back())) {
          // append the key of p
          pair_after_filtering.push_back(p.first);
        }

        // we only want up to MAX_SITE_PER_SLICE combinations
        if (pair_after_filtering.size() >= this->parameters.MAX_SITE_PER_SPLICE) {
          // so just disreagard any after that point
          break;
        }
      }
    }

    int support_count_total = 0;

    for (const auto & p : pair_after_filtering) {
      // now we want to add the values of keys matching certain criteria
      for (const auto & it : tmp_pair) {
        if ((abs(it.first.start - p.start) + abs(it.first.end - p.end)) < this->parameters.MAX_TS_DIST) {
          support_count_total += it.second;
        }
      }
    }

    if ((float(support_count_total) / this->lr_pair.size()) <= this->parameters.MIN_SUP_PCT) {
      // not enough support counts - this often happens in the case when the TSS/TES is strongly degraded

      if (pair_enrich.size() == 0) {
        for (const auto & p : tmp_pair) {
          // add it to pair_enrich
          // we need to sum up all of the elements in tmp_pair that meet our criteria
          int sum = 0;
          for (const auto & it : tmp_pair) {
            if ((abs(it.first.start - p.first.start) + abs(it.first.end - p.first.end)) < this->parameters.MAX_TS_DIST) {
              sum += it.second;
            }
          }
          pair_enrich.push_back({p.first, sum});
        }
        std::sort(
          pair_enrich.begin(), pair_enrich.end(),
          [] (const auto & p1, const auto & p2) {
            // sign is flipped because we are reverse sorting
            return (p1.second < p2.second);
          }
        );
      }

      std::vector<int>
      tmp_ex;
      tmp_ex.push_back(int(pair_enrich[0].first.start));
      for (int i : pair.first) {
        tmp_ex.push_back(i);
      }
      tmp_ex.push_back(int(pair_enrich[0].first.end));

      this->raw_isoforms[tmp_ex] = pair.second.size();
      this->strand_counts[tmp_ex] = this->strand_counts[pair.first];
    } else {
      // now, add filtered TSS/TES to raw_isoforms
      for (const auto & p : pair_after_filtering) {
        std::vector<int>
        tmp_ex;
        tmp_ex.push_back(int(p.start));
        for (int i : pair.first) {
          tmp_ex.push_back(i);
        }
        tmp_ex.push_back(int(p.end));

        // add up the values to assign to raw_isoforms
        int sum = 0;
        for (const auto & it : tmp_pair) {
          if (abs(it.first.start - p.start) + abs(it.first.end - p.end) < this->parameters.MAX_TS_DIST) {
            sum += it.second;
          }
        }

        this->raw_isoforms[tmp_ex] = sum;
        this->strand_counts[tmp_ex] = this->strand_counts[pair.first];
      }
    }
  }
}


/*
  this function does most of the heavy lifting in the Isoform class
  it does of stuff - honestly my eyes glaze over when I try to understand it all at once
*/
void Isoforms::match_known_annotation 
(
  std::unordered_map<std::string, Junctions> transcript_to_junctions,
  std::unordered_map<std::string, Pos> transcript_dict,
  std::unordered_map<std::string, std::vector<StartEndPair>> gene_dict,
  GeneBlocks one_block,
  std::unordered_map<std::string, std::string> fa_dict
)
{

  std::set<int>
  splice_site = get_splice_site(transcript_to_junctions, one_block.transcript_list);

  std::vector<std::vector<int>>
  junction_list;
  std::unordered_map<std::vector<int>, std::string>
  junction_dict;

  std::vector<std::vector<int>>
  exons_list;
  std::unordered_map<std::vector<int>, std::string>
  exons_dict;

  // populate junction_list, junction_dictionary
  for (const auto & i : one_block.transcript_list) {
    junction_list.push_back(transcript_to_junctions[i].junctions);
    junction_dict[transcript_to_junctions[i].junctions] = i;

    // add the new entry to exons_list
    exons_list.push_back({transcript_to_junctions[i].left[0]});
    exons_list.back().insert(
      exons_list.back().end(),
      transcript_to_junctions[i].junctions.begin(), 
      transcript_to_junctions[i].junctions.end()
    );
    exons_list.back().push_back(transcript_to_junctions[i].right[0]);

    exons_dict[exons_list.back()] = i;
  }

  Junctions
  TSS_TES_site = get_TSS_TES_site(transcript_to_junctions, one_block.transcript_list);

  // if this is a single-exon read, we're done
  if (splice_site.size() == 0) {
    return;
  }

  for (const auto & [exon_key, exon_val] : this->single_block_dict) {
    for (const auto & i : one_block.transcript_list) {
      int tmp_std;

      if (this->parameters.STRAND_SPECIFIC == 0) {
        // we need to convert the pair to a vector for this lookup
        tmp_std = this->strand_counts[{exon_key.start, exon_key.end}]; // this is where i'm up to
      } else {
        // little shorthand if else for you there ;)
        transcript_dict[i].strand == '+' ? tmp_std = 1 : tmp_std = 0;
      }

      if ((transcript_to_junctions[i].junctions.size() == 0) &&
          (tmp_std == this->strand_counts[{exon_key.start, exon_key.end}])) {
        if ((abs(exon_key.start - transcript_to_junctions[i].left[0]) < this->parameters.MAX_TS_DIST) &&
            (abs(exon_key.end - transcript_to_junctions[i].right[0]) < this->parameters.MAX_TS_DIST)) {
          std::pair<int, int>
          known_exons = {transcript_to_junctions[i].left[0], transcript_to_junctions[i].right[0]};
          if (this->known_isoforms.count({known_exons.first, known_exons.second})) { // check if it is already known
            // it's already known, just update the entry
            this->known_isoforms[std::vector<int>{known_exons.first, known_exons.second}] = {
              long(this->known_isoforms[{known_exons.first, known_exons.second}].support_count + this->single_blocks[this->single_block_dict[exon_key]].size()),
              i,
              transcript_dict[i].parent_id
            };
          } else {
            // it's totally new, create a fresh entry for it
            this->known_isoforms[{known_exons.first, known_exons.second}] = {
              long(this->single_blocks[this->single_block_dict[exon_key]].size()),
              i,
              transcript_dict[i].parent_id
            };
            if (!this->parameters.STRAND_SPECIFIC) {
              this->strand_counts[{known_exons.first, known_exons.second}] = 1;
              if (transcript_dict[i].strand == '-') {
                this->strand_counts[{known_exons.first, known_exons.second}] = -1;
              }
            }
          }
        }
      }
    }
  }

  for (auto const & [raw_iso_key, raw_iso_val] : this->raw_isoforms) {
    // line 730 in sc_longread.py
    int found = 0;

    // we need to check if there are any repeated values in raw_iso_key
    std::set<int>
    raw_iso_key_set;
    for (const auto & i : raw_iso_key) {
      raw_iso_key_set.insert(i);
    }

    // ignore raw_iso data that has repeated entries
    if (raw_iso_key.size() > raw_iso_key_set.size()) {
      continue;
    }

    int tmp_std;

    for (const auto & i : one_block.transcript_list) {
      if (!this->parameters.STRAND_SPECIFIC) {
        tmp_std = this->strand_counts[raw_iso_key];
      } else {
        transcript_dict[i].strand == '+' ? tmp_std = 1 : tmp_std = -1;
      }

      // same number of exons, same strand
      if ((raw_iso_key.size() - 2 == transcript_to_junctions[i].junctions.size()) &&
          (tmp_std == this->strand_counts[raw_iso_key])) {
        int iso_is_within_max_dist = 1;
        for (int j = 0; j < std::min(raw_iso_key.size() - 1, transcript_to_junctions[i].junctions.size()); j++) {
          if (abs(raw_iso_key[j+1] - transcript_to_junctions[i].junctions[j]) > this->parameters.MAX_DIST) {
            // we don't want any to be greater than MAX_DIST
            iso_is_within_max_dist = 0;
            break;
          }
        }

        if (iso_is_within_max_dist) {
          if ((abs(raw_iso_key.front() - transcript_to_junctions[i].left[0]) < this->parameters.MAX_DIST) &&
              (abs(raw_iso_key.back() - transcript_to_junctions[i].right[0]) < this->parameters.MAX_DIST)) {
            // populate known_exons with transcript_to_junctions
            std::vector<int>
            known_exons = {transcript_to_junctions[i].left[0]};
            for (const auto & j : transcript_to_junctions[i].junctions) {
              known_exons.push_back(j);
            }
            known_exons.push_back(transcript_to_junctions[i].right[0]);

            if (this->known_isoforms.count(known_exons)) {
              // line 747 sc_longread.py
              
              this->known_isoforms[known_exons] = {
                long(std::max(this->known_isoforms[known_exons].support_count, long(raw_iso_val))),
                i,
                transcript_dict[i].parent_id
              };
            } else {
              this->known_isoforms[known_exons] = {
                raw_iso_val,
                i,
                transcript_dict[i].parent_id
              };

              if (!this->parameters.STRAND_SPECIFIC) {
                // if it's not strand specific protocal, use annotation
                this->strand_counts[known_exons] = 1;
                if (transcript_dict[i].strand == '-') {
                  this->strand_counts[known_exons] = -1;
                } else {
                  this->strand_counts[known_exons] = this->strand_counts[raw_iso_key];
                }
              }

              found = 1;
              break;
            }
          }
        }
      }
    }

    if (!found) {
      std::vector<int>
      new_exons = find_best_splice_chain(raw_iso_key, junction_list, this->parameters.MAX_SPLICE_MATCH_DIST);

      int is_strictly_increasing = 1;
      for (int i = 0; i < new_exons.size() - 1; i++) {
        if (!(new_exons[i] < new_exons[i + 1])) {
          // then the function is not strictly increasing
          is_strictly_increasing = 0;
          break;
        }
      }

      if (!is_strictly_increasing) {
        new_exons = raw_iso_key;
      }

      for (int i = 0; i < new_exons.size(); i++) {
        int a_site = new_exons[i];

        if (i == 0) { // the left-most site
          int closest = take_closest(TSS_TES_site.left, a_site);

          if ((abs(closest - a_site) < this->parameters.MAX_TS_DIST) && 
              (abs(closest < raw_iso_key[i + 1]))) {
            new_exons[i] = closest;
          } else {
            new_exons[i] = a_site;
          }
        } else if (0 < i < raw_iso_key.size() - 1) { // a site from the middle somewhere
          std::vector<int>
          splice_site_vec;
          for (const auto & site : splice_site) {
            splice_site_vec.push_back(site);
          }
          int closest = take_closest(splice_site_vec, a_site);

          if ((abs(closest - a_site) < this->parameters.MAX_SPLICE_MATCH_DIST) &&
              (closest > new_exons.back()) &&
              (closest < raw_iso_key[i + 1])) {
            new_exons[i] = closest;
          } else {
            new_exons[i] = a_site;
          }
        } else { // then it must be the right-most site
          int closest = take_closest(TSS_TES_site.right, a_site);
          if ((abs(closest - a_site) < this->parameters.MAX_TS_DIST) &&
              (closest > new_exons.back())) {
            new_exons[i] = closest;
          } else {
            new_exons[i] = a_site;
          }
        }
      }


      is_strictly_increasing = 1;
      for (int i = 0; i < new_exons.size() - 1; i++) {
        if (!(new_exons[i] < new_exons[i - 1])) {
          // then the vector is not strictly increasing
          is_strictly_increasing = 0;
          break;
        }
      }

      if (is_strictly_increasing) {
        if (this->new_isoforms.count(new_exons) == 0) {
          this->new_isoforms[new_exons] = {this->raw_isoforms[raw_iso_key], "", ""};
          this->strand_counts[new_exons] = this->strand_counts[raw_iso_key];
        } else {
          this->new_isoforms[new_exons] = {this->new_isoforms[new_exons].support_count + raw_iso_val, "", ""};
        }
      }
    }
  }

  // line 800 of sc_longread.py

  // remove incomplete transcript (due to 3' bias)
  if (this->parameters.REMOVE_INCOMP_READS) {
    std::vector<std::vector<int>>
    delete_key;

    for (const auto & [new_isoform_key, new_isoform_val] : this->new_isoforms) {
      // 1. match to known isoform detected and remove new isoform if within known isoform
      if (this->known_isoforms.count(new_isoform_key)) {
        delete_key.push_back(new_isoform_key);
        continue;
      }

      for (const auto & [known_isoform_key, known_isoform_val] : this->known_isoforms) {
        if ((new_isoform_key.size() < known_isoform_key.size()) &&
            (if_exon_contains(known_isoform_key, std::vector<int>(new_isoform_key.begin()+2, new_isoform_key.end()), 1))) {
          std::vector<char>
          s_l = std::vector<char>( // get a slice of 15 ints just before new_isoform_key[0]
            fa_dict[this->ch].begin() + (new_isoform_key[0] - 15), 
            fa_dict[this->ch].begin() + (new_isoform_key[0])
          );

          std::vector<char>
          s_r = std::vector<char>( // get a slice of 15 ints just after new_isoform_key[0]
            fa_dict[this->ch].begin() + new_isoform_key[0],
            fa_dict[this->ch].begin() + new_isoform_key[0] + 15
          );

          if ((std::count(s_l.begin(), s_l.end(), 'T') > 10) ||
              (std::count(s_l.begin(), s_l.end(), 'A') > 10) ||
              (std::count(s_r.begin(), s_r.end(), 'T') > 10) ||
              (std::count(s_r.begin(), s_r.end(), 'A') > 10)) {
            // remove keys with too many Ts or As at the beginning or end
            delete_key.push_back(new_isoform_key);
            break;
          }
        }
        else if ((new_isoform_key.size() < known_isoform_key.size()) &&
                (if_exon_contains(known_isoform_key, std::vector<int>(new_isoform_key.begin(), new_isoform_key.end() - 2), 1))) {
          std::vector<char>
          s_l = std::vector<char>( // get a slice of 15 ints just before new_isoform_key.back()
            fa_dict[this->ch].begin() + (new_isoform_key.back() - 15), 
            fa_dict[this->ch].begin() + (new_isoform_key.back())
          );

          std::vector<char>
          s_r = std::vector<char>( // get a slice of 15 ints just after new_isoform_key.back()
            fa_dict[this->ch].begin() + new_isoform_key.back(),
            fa_dict[this->ch].begin() + new_isoform_key.back() + 15
          );
          
          if ((std::count(s_l.begin(), s_l.end(), 'T') > 10) ||
              (std::count(s_l.begin(), s_l.end(), 'A') > 10) ||
              (std::count(s_r.begin(), s_r.end(), 'T') > 10) ||
              (std::count(s_r.begin(), s_r.end(), 'A') > 10)) {
            // remove keys with too many Ts or As at the beginning or end
            delete_key.push_back(new_isoform_key);
            break;
          }
        }

        if ((new_isoform_key.size() < known_isoform_key.size()) && if_exon_contains(known_isoform_key, new_isoform_key, this->parameters.MAX_TS_DIST)) {
          float sim_pct_sq = pow(get_exon_sim_pct(new_isoform_key, known_isoform_key), 2);
          if (new_isoform_val.support_count < (1 + sim_pct_sq * this->parameters.REMOVE_INCOMP_READS) * this->known_isoforms[known_isoform_key].support_count) {
            // remove keys with too much similarity in coverage
            delete_key.push_back(new_isoform_key);
            break;
          }
        }

        // add to delete keys if they are the same
        if (new_isoform_key == known_isoform_key) {
          delete_key.push_back(new_isoform_key);
          break;
        }
      }
    }

    // convert delete_key to a set
    std::set<std::vector<int>>
    delete_key_set;
    for (const auto & v : delete_key) {
      delete_key_set.insert(v);
    }

    // entirely remove the keys in delete_key_set from new_isoforms 
    for (const auto & key : delete_key_set) {
      this->new_isoforms.erase(key);
    }

    // reset delete key
    delete_key = {};

    // 2. match to known isoform in annotation, and remove new_isoform if it's very similar to annotation
    for (const auto & [new_isoform_key, new_isoform_val] : this->new_isoforms) {
      for (const auto & [known_isoform_key, known_isoform_va] : exons_dict) {
        if (!this->known_isoforms.count(known_isoform_key)) {
          if (new_isoform_key.size() < known_isoform_key.size() &&
              if_exon_contains(known_isoform_key, new_isoform_key, this->parameters.MAX_TS_DIST)) {
            auto 
            sim_pct = get_exon_sim_pct(new_isoform_key, known_isoform_key);

            if (sim_pct > 0.95) {
              delete_key.push_back(new_isoform_key);
              this->known_isoforms[known_isoform_key] = {
                this->new_isoforms[new_isoform_key].support_count,
                exons_dict[known_isoform_key],
                transcript_dict[exons_dict[known_isoform_key]].parent_id
              };

              this->strand_counts[known_isoform_key] = 1;
              if (transcript_dict[exons_dict[known_isoform_key]].strand != '+') {
                this->strand_counts[known_isoform_key] = -1;
              }
            }
          } else if (new_isoform_key == known_isoform_key) {
            delete_key.push_back(new_isoform_key);
            this->known_isoforms[known_isoform_key] = {
              this->new_isoforms[new_isoform_key].support_count,
              exons_dict[known_isoform_key],
              transcript_dict[exons_dict[known_isoform_key]].parent_id
            };

            this->strand_counts[known_isoform_key] = 1;
            if (transcript_dict[exons_dict[known_isoform_key]].strand != '+') {
              this->strand_counts[known_isoform_key] = -1;
            }
          }
        }
      }
    }

    // reset delete key set
    delete_key_set = {};
    for (const auto & v : delete_key) {
      delete_key_set.insert(v);
    }

    // entirely remove the keys in delete_key_set from new_isoforms 
    for (const auto & key : delete_key_set) {
      this->new_isoforms.erase(key);
    }

    // 3. Match among new isoform
    if (this->new_isoforms.size() > 1) {
      // reset delete key
      delete_key = {};

      // make a vector vector isos, and populate it with isoform keys
      std::vector<std::vector<int>> 
      isos;
      for (const auto & [i, j] : this->new_isoforms) {
        isos.push_back(i);
      }

      for (int i = 0; i<isos.size() - 1; i++) {
        for (int j = i+1; j < isos.size(); j++) {
          if (isos[i].size() < isos[j].size() &&
              if_exon_contains(isos[j], isos[i], this->parameters.MAX_TS_DIST)) {
            auto 
            sim_pct_sq = pow(get_exon_sim_pct(isos[j], isos[i]), 2);
            
            if (this->new_isoforms[isos[i]].support_count < sim_pct_sq * this->parameters.REMOVE_INCOMP_READS * this->new_isoforms[isos[j]].support_count) {
              delete_key.push_back(isos[i]);
            }
          } else if (isos[i].size() > isos[j].size() && 
                      if_exon_contains(isos[i], isos[j], this->parameters.MAX_TS_DIST)) {
            auto 
            sim_pct_sq = pow(get_exon_sim_pct(isos[j], isos[i]), 2);

            if (this->new_isoforms[isos[i]].support_count < sim_pct_sq * this->parameters.REMOVE_INCOMP_READS * this->new_isoforms[isos[j]].support_count) {
              delete_key.push_back(isos[j]);
            }
          }
        }
      }

      // convert delete_key to a set
      std::set<std::vector<int>>
      delete_key_set;
      for (const auto & v : delete_key) {
        delete_key_set.insert(v);
      }

      // entirely remove the keys in delete_key_set from new_isoforms 
      for (const auto & key : delete_key_set) {
        this->new_isoforms.erase(key);
      }
    }
  }


  // match to gene

  std::vector<std::vector<int>>
  delete_key;
  std::map<std::vector<int>, Iso>
  update_iso_dict;

  if (this->new_isoforms.size() > 0) {
    for (const auto & [isoform_key, isoform_val] : this->new_isoforms) {
      // calculate a sum of the pairs of the isoform
      int iso_len = 0;
      auto iso_pairs = pairwise(isoform_key);
      for (const auto & i : iso_pairs) {
        iso_len += i.end - i.start;
      }
      
      // make tmp to store exon overlap and gene name, and populate it
      std::vector<std::pair<int, std::string>>
      tmp;
      for (const auto & [ge, tr] : one_block.gene_to_tr) {
        tmp.push_back({exon_overlap(isoform_key, gene_dict[ge]), ge});
      }

      if (this->parameters.STRAND_SPECIFIC) { // if it has strand-specific protocol, use read
        char stnd = this->strand_counts[isoform_key] == 1 ? '+' : '-';

        // update tmp, only including values which match the sign
        std::vector<std::pair<int, std::string>>
        new_tmp;
        for(const auto & it : tmp) {
          if (transcript_dict[one_block.gene_to_tr[it.second][0]].strand == stnd) {new_tmp.push_back(it);}
        }
        tmp = new_tmp;

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
        if (isoform_key.front() >= gene_dict[tmp[0].second].front().start &&
            isoform_key.back() <= gene_dict[tmp[0].second].back().end) {
          this->new_isoforms[isoform_key] = {
            this->new_isoforms[isoform_key].support_count,
            "",
            tmp[0].second
          };
        }
        else if ((exon_overlap(std::vector<int>(isoform_key.begin() + 2, isoform_key.begin() + 4), gene_dict[tmp[0].second]) > 0) &&
                 (exon_overlap(std::vector<int>(isoform_key.end() - 4, isoform_key.end() - 2), gene_dict[tmp[0].second]) > 0)) {
          if ((isoform_key[1]-isoform_key[0] <= this->parameters.MIN_FL_EXON_LEN) && 
              (exon_overlap(std::vector<int>(isoform_key.begin() + 2, isoform_key.begin() + 4), gene_dict[tmp[0].second]) == 0)) {
            auto
            ba = std::vector<char>(fa_dict[this->ch].begin() + isoform_key[0], fa_dict[this->ch].begin() + isoform_key[1]);
            
            if (std::count(ba.begin(), ba.end(), 'T') > 0.7 * ba.size() || std::count(ba.begin(), ba.end(), 'A') > 0.7 * ba.size()) {  // leftover polyA tail 
              delete_key.push_back(isoform_key);
              update_iso_dict[std::vector<int>(isoform_key.begin()+2, isoform_key.end())] = {
                this->new_isoforms[isoform_key].support_count,
                "",
                tmp[0].second
              };
            }
          } else if ((isoform_key.back() - isoform_key.rbegin()[1] <= this->parameters.MIN_FL_EXON_LEN) &&
                    (exon_overlap(std::vector<int>(isoform_key.begin(), isoform_key.begin() + 2), gene_dict[tmp[0].second])== 0)) {
            auto
            ba = std::vector<char>(fa_dict[this->ch].begin() + isoform_key.rbegin()[1], fa_dict[this->ch].begin() + isoform_key.back());

            if (std::count(ba.begin(), ba.end(), 'T') > 0.7 * ba.size() || std::count(ba.begin(), ba.end(), 'A') > 0.7 * ba.size()) {
              delete_key.push_back(isoform_key);
              update_iso_dict[std::vector<int>(isoform_key.begin(), isoform_key.end() - 1)] = {
                this->new_isoforms[isoform_key].support_count,
                "",
                tmp[0].second
              };
            }
          } else {
            this->new_isoforms[isoform_key] = {this->new_isoforms[isoform_key].support_count, "", tmp[0].second};
          }
        } else {
          if (tmp[0].first > 0.8 * iso_len) {
            if ((isoform_key[1] - isoform_key[0] < this->parameters.MIN_FL_EXON_LEN) || 
                (isoform_key.rbegin()[0] - isoform_key.rbegin()[1] < this->parameters.MIN_FL_EXON_LEN)) {
              continue; // alignment artifact
            } else {
              // might be real eRNA
              this->new_isoforms[isoform_key] = {
                this->new_isoforms[isoform_key].support_count,
                "",
                tmp[0].second
              };
            }
          } else {
            // get the total of elements in isoform_key that are also in splice_count
            int total = 0; 
            for (const auto & it : isoform_key) {
              if (splice_site.count(it) > 0) {
                total++;
              }
            }

            if (total > 4) {
              // more than 4 splice site match
              this->new_isoforms[isoform_key] = {
                this->new_isoforms[isoform_key].support_count,
                "",
                tmp[0].second
              };
            }
          }
        }
      }
    }

    // line 910
    // remove all the delete_key entries from new_isoforms
    if (delete_key.size() > 0) {
      for (const auto & key : delete_key) {
        this->new_isoforms.erase(key);
      }
    }

    for (const auto & [iso_key, iso_val] : this->new_isoforms) {
      // skip if there are no matching genes
      if (iso_val.gene_id == "") {
        continue;
      }

      // add iso_val.gene_id to ge_dict if it's not already there
      if (this->ge_dict.count(iso_val.gene_id) == 0) {
        this->ge_dict[iso_val.gene_id] = {};
      }
    
      // and then update the entry with iso_key
      this->ge_dict[iso_val.gene_id].push_back(iso_key);

      if (!this->parameters.STRAND_SPECIFIC) {
        if (this->strand_counts[iso_key] == 0) {
          this->strand_counts[iso_key] = 1;
          if (transcript_dict[one_block.gene_to_tr[this->new_isoforms[iso_key].gene_id][0]].strand == '-') {
            this->strand_counts[iso_key] = -1;
          }
        }
      }
    }
  }

  for (const auto & [iso_key, iso_val] : this->known_isoforms) {
    // add the iso_val's gene_id to ge_dict if it's not already there
    if (this->ge_dict.count(iso_val.gene_id) == 0) {
      this->ge_dict[iso_val.gene_id] = {};
    }
    // then update the value with the new iso_key
    this->ge_dict[this->known_isoforms[iso_key].gene_id].push_back(iso_key);
  }
}

std::string
Isoforms::isoform_to_gtt3(int isoform_pct=-1) {
  std::vector<std::string>
  gff_rec = {};

  std::map<std::string, int>
  transcript_id_dict;

  if (this->new_isoforms.size() + this->known_isoforms.size() == 0) {
    return "";
  }

  for (const auto & [g_key, g_val] : this->ge_dict) {
    std::vector<std::string>
    gff_tmp = {};

    // add up all the support_count from all the genes in g_val
    int total_cnt = 0;
    for (const auto & e : g_val) {
      if (this->new_isoforms.count(e) > 0) {total_cnt += this->new_isoforms[e].support_count;}
      if (this->known_isoforms.count(e) > 0) {total_cnt += this->known_isoforms[e].support_count;}
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
    gff_entry << this->ch << '\t' // _ch
              << "." << '\t' // _sr
              << "gene" << '\t' // _ty
              << *std::min_element(firsts.begin(), firsts.end()) + 1 << '\t' // _st
              << *std::max_element(lasts.begin(), lasts.end()) << '\t' // _en
              << "." << '\t' // _sc
              << (this->strand_counts[g_val[0]] == 1 ? '+' : '-') << '\t' // _stnd
              << "." << "\t" // _ph
              << "ID=gene:" << g_key << ";gene_id=" << g_key << ";support_count=" << total_cnt; // _attr
    
    gff_tmp.push_back(gff_entry.str());

    for (const auto & exons : g_val) {
      if (this->new_isoforms.count(exons) && this->known_isoforms.count(exons)) {
        std::cout << "BOTH in new and known\n";
      }

      std::string 
      source;
      long 
      support_count;
      std::stringstream 
      tp_id;

      if (this->new_isoforms.count(exons)) {
        source = "new";
        support_count = this->new_isoforms[exons].support_count;

        if (0 < isoform_pct && isoform_pct < 1 && support_count < (isoform_pct * total_cnt)) {
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

      if (support_count < this->parameters.STRAND_SPECIFIC) {
        continue;
      }
      auto exon_idx = 1;

      std::stringstream gff_entry;
      gff_entry << this->ch << '\t' // _ch // chromosome
                << source << '\t' // _sr // source
                << "transcript" << '\t' // _ty // 
                << exons[0] + 1 << '\t' // _st // start
                << exons.back() << '\t' // _en // end
                << "." << '\t' // _sc 
                << (this->strand_counts[exons] == 1 ? '+' : '-') << '\t' // _stnd // strand
                << "." << "\t" // _ph
                << "ID=transcript:" << tp_id.str() << ";transcript_id=" << tp_id.str() << ";Parent=gene:" << g_key << ";support_count=" << support_count << ";source=" << source; // _attr
      gff_tmp.push_back(gff_entry.str());
      exon_idx++;
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
    output << "\n";
    return output.str();
  
  }

  return "";
}

