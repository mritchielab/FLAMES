#include <map>
#include <vector>
#include <algorithm>
#include <cmath>
#include <numeric>

#include <iostream>
#include <fstream>

#include "junctions.hpp"
#include "GeneBlocks.hpp"
#include "misc.hpp"

class Isoforms
{
  private:
    // these values will all be extracted from config
    const int MAX_DIST, MAX_TS_DIST, MAX_SPLICE_MATCH_DIST,
              MAX_SITE_PER_SLICE, MIN_SUP_CNT, MIN_FL_EXON_LEN,
              MIN_SUP_PCT, STRAND_SPECIFIC, REMOVE_INCOMP_READS;
  public:
    std::string ch;
    
    std::map<std::vector<int>, int> 
    junction_dict;
    std::vector<std::vector<std::vector<int>>> 
    junction_list;

    std::map<std::vector<int>, std::vector<std::pair<int, int>>>
    lr_pair;
    std::vector<int> 
    left;
    std::vector<int> 
    right;

    std::map<std::pair<int, int>, int> 
    single_block_dict;
    std::vector<std::vector<std::pair<int, int>>> 
    single_blocks;
    std::map<std::vector<int>, std::vector<int>> 
    strand_cnt;
    std::map<std::string, int> 
    new_isoforms;
    std::map<std::vector<int>, Iso> 
    known_isoforms;
    std::map<std::vector<int>, std::vector<std::pair<int, int>>> 
    raw_isoforms;
    std::map<std::string, int> 
    ge_dict;

    Isoforms(std::string ch, std::map<std::string, int> config);
    void add_isoform(Junctions junctions, bool is_reversed);
    void add_one(Junctions junctions, bool strand);
    void update_one(Junctions junctions, std::vector<int> key, bool strand);
    int len();
    void update_all_splice();

    void filter_TSS_TES(std::ofstream out_f, Junctions known_site, float fdr_cutoff); 

    //unused
    std::pair<std::vector<int>, std::map<int, int>>
    group_sites(std::vector<int> l, int smooth_window, int min_threshold);

    void Isoforms::match_known_annotation (
  std::map<std::string, Junctions> transcript_to_junctions,
  std::map<std::string, Pos> transcript_dict,
  std::map<std::string, int> gene_dict,
  GeneBlocks one_block,
  std::map<std::string, int> fa_dict
  );

};

Isoforms::Isoforms(std::string ch, std::map<std::string, int> config)
: MAX_DIST              (config["MAX_DIST"]),
  MAX_TS_DIST           (config["MAX_TS_DIST"]),
  MAX_SPLICE_MATCH_DIST (config["MAX_SPLICE_MATCH_DIST"]),
  MAX_SITE_PER_SLICE    (config["Max_site_per_splice"]),
  MIN_SUP_CNT           (config["Min_sup_cnt"]),
  MIN_FL_EXON_LEN       (config["min_fl_exon_len"]),
  MIN_SUP_PCT           (config["Min_sup_pct"]),
  STRAND_SPECIFIC       (config["strand_specific"]),
  REMOVE_INCOMP_READS   (config["remove_incomp_reads"])
{
  /*
    initialises the object,
    extracting all the const values from config
    
    (why are config keys a mix of CAPS, lower, and Sentence Case?)
    (I have no idea, that's how they were when I found them)
    (change it later)
  */
  
  this->ch = ch;
}

void Isoforms::add_isoform(Junctions junctions, bool is_reversed)
{
  if (junctions.junctions.size() == 0) // single exon reads
  {
    if (this->single_block_dict.size() == 0)
    {
      // this will be the first thing in single_block_dict
      // so, nothing to check
      this->add_one(junctions, is_reversed);
    }
    else
    {
      // there is already stuff in single_block_dict
      // first, check if this key is already in there
      std::pair<int, int> target_key = 
      {junctions.left[0], junctions.right[0]};

      if (this->single_block_dict.count(target_key) > 0)
      {
        // this time the key must be converted to vector
        // we really need to fix this
        std::vector<int> target_key_vector =
        {target_key.first, target_key.second};

        // the key is already present - so just update it
        this->update_one(junctions, target_key_vector, is_reversed);
      }
      else
      {
        // the key is not present in the dict yet
        // check if there's anything very close to the key
        bool found = false;

        for (auto const& [current_key, current_value] : this->single_block_dict)
        {
          if ((abs(current_key.first - target_key.first) < this->MAX_TS_DIST) and
              (abs(current_key.second - target_key.second) < this->MAX_TS_DIST))
          {
            // we've found one that's very similar
            // just update it
            this->update_one(junctions, {current_key.first, current_key.second}, is_reversed);
            found = true;
            break;
          }
        }

        if (!found)
        {
          // we didn't find any similar ones
          // so just add it
          this->add_one(junctions, is_reversed);
        }
      }
    }
  }
  else // multi-exon reads
  {
    if (this->junction_dict.size() == 0)
    {
      // there's nothing in the junction dict yet
      // just add it
      this->add_one(junctions, is_reversed);
    }
    else
    {
      // there's stuff in junction_dict
      // check to see if there's an exact match

      if (this->junction_dict.count(junctions.junctions) == 0)
      {
        // the key is not present in the dict yet
        // just add it
        this->add_one(junctions, is_reversed);
      }
      else
      {
        // the key is not directly present in the dict
        // check for similar keys

        bool found = false;

        for (auto & [current_key, current_value] : this->junction_dict)
        {
          if (junctions.junctions.size() == current_key.size())
          {
            // they are the same size. see if they are similar
            
            bool similar = true;
            for (int i = 0; i < current_key.size(); i++)
            {
              if (abs(junctions.junctions[i] - current_key[i]) > this->MAX_DIST)
              {
                // they are too far apart
                similar = false;
                break;
              }
            }

            if (similar)
            {
              // if they are still similar,
              // then we have found an entry that is close enough
              
              this->update_one(junctions, current_key, is_reversed);
              found = true;
              break;
            }
          }
        }
        
        if (!found)
        {
          // we went through and found nothing similar enough
          // so just add it to the dict
          this->add_one(junctions, is_reversed);
        }
      }
    }
  }
};

void Isoforms::add_one(Junctions junctions, bool strand)
{
  if (junctions.junctions.size() == 0) // single-exon
  {
    std::pair<int, int> key = 
    {junctions.left[0], junctions.right[0]};

    this->single_block_dict[key] = this->single_blocks.size();
    this->single_blocks.push_back({key});
    this->strand_cnt[std::vector<int> {key.first, key.second}] = {};

    int multiplier = 1;
    if (strand) {multiplier = -1;}
    this->strand_cnt[std::vector<int> {key.first, key.second}].push_back(multiplier * this->STRAND_SPECIFIC);
  }
  else // multi-exon reads
  {
    this->junction_dict[junctions.junctions] = this->junction_list.size();
    this->junction_list.push_back({junctions.junctions});
    this->left.push_back(junctions.left[0]);
    this->right.push_back(junctions.right[0]);
    this->lr_pair[junctions.junctions] = {};
    
    this->strand_cnt[junctions.junctions] = {};
    int multiplier = 1;
    if (strand) {multiplier = -1;}
    this->strand_cnt[junctions.junctions].push_back(multiplier * this->STRAND_SPECIFIC);
  }
}

void Isoforms::update_one(Junctions junctions, std::vector<int> key, bool strand)
{
  if (junctions.junctions.size()==0) // single-exon reads
  {
    std::pair<int, int> new_element = 
    {junctions.left[0], junctions.right[0]};

    // we need to express the key as a pair here,
    // because we're looking through single_block_dict which is full of pairs
    // honestly this is gonna become a big problem with this 
    std::pair<int, int> key_pair =
    {key[0], key[1]};

    this->single_blocks[this->single_block_dict[key_pair]].push_back(new_element);
  }
  else // multi-exon reads
  {
    this->junction_list[this->junction_dict[key]].push_back(junctions.junctions);
    this->left.push_back(junctions.left[0]);
    this->right.push_back(junctions.right[0]);
    this->lr_pair[key].push_back(std::pair<int, int> {junctions.left[0], junctions.right[0]});
  }
  int multiplier = 1;
  if (this->STRAND_SPECIFIC) {multiplier = -1;}  
  this->strand_cnt[key].push_back(multiplier * this->STRAND_SPECIFIC);
}

int Isoforms::len()
{
  return (this->junction_dict.size() + this->single_block_dict.size());
}

void Isoforms::update_all_splice()
{

  std::map<std::vector<int>, int>
  junction_tmp;
  std::map<std::pair<int, int>, int>
  single_block_tmp;
  std::map<std::vector<int>, std::vector<int>>
  strand_cnt_tmp;
  std::map<std::vector<int>, std::vector<std::pair<int,int>>>
  lr_pair_tmp;

  // look through junction_dict
  for (auto const& [key, val] : this->junction_dict)
  {
    if (this->junction_list[val].size() >= this->MIN_SUP_CNT)
    {
      // this will store the most common junctions
      // variable to store the counts
      std::map<std::vector<int>, int>
      junction_list_counts;

      // populate the counts container
      for (auto i : this->junction_list[val])
      {
        junction_list_counts[i] ++;
      }
      
      // pull out the key for the 
      // most common value in the counts map
      auto new_key = (*std::max_element
      (
        junction_list_counts.cbegin(), junction_list_counts.cend()
      )).first;

      junction_tmp[new_key] = this->junction_dict[key];
      
      std::map<int, int>
      strand_cnt_counts;
      for (auto i : this->strand_cnt[key])
      {
        strand_cnt_counts[i] ++;
      }
      // get the key for the most common
      // value in max_strand_cnt
      auto max_strand_cnt = (*std::max_element
      (
        strand_cnt_counts.cbegin(), strand_cnt_counts.cend()
      )).first;

      strand_cnt_tmp[new_key] = {max_strand_cnt};
      lr_pair_tmp[new_key] = this->lr_pair[key];
    }
  }

  // look through single_block_dict
  for (auto const& [key, val] : this->single_block_dict)
  {
    if (this->single_blocks[val].size() >= this->MIN_SUP_CNT)
    {
      // this will store the counts
      std::map<std::pair<int, int>, int>
      single_blocks_counts;
      for (auto i : this->single_blocks[val])
      {
        single_blocks_counts[i] ++;
      }

      auto new_key = (*std::max_element
      (
        single_blocks_counts.cbegin(), single_blocks_counts.cend()
      )).first;

      single_block_tmp[new_key] = this->single_block_dict[key];
      std::map<int, int>
      strand_cnt_counts;
      for (auto i : this->strand_cnt[{key.first, key.second}])
      {
        strand_cnt_counts[i] ++;
      }
      auto max_strand_cnt = (*std::max_element(
        strand_cnt_counts.cbegin(), strand_cnt_counts.cend()
      )).first;
      strand_cnt_tmp[{new_key.first, new_key.second}] = {max_strand_cnt};
    }
  }

  this->single_block_dict = single_block_tmp;
  this->junction_dict = junction_tmp;
  this->strand_cnt = strand_cnt_tmp;
  this->lr_pair = lr_pair_tmp;
}

void Isoforms::filter_TSS_TES(std::ofstream out_f, Junctions known_site={}, float fdr_cutoff=0.01)
{
  std::string bedgraph_fmt = "{_ch}\t{_st}\t{_en}\t{_sc}\n";

  auto filter_site = [fdr_cutoff] (std::map<int, int> list_counts)
  {
    // to sort the map, we need it as a vector first.
    // note: i should find a better method for doing this

    std::vector<std::pair<int, int>>
    mx;

    for (const auto & i : list_counts)
    {
      mx.emplace_back(i);
    }

    std::sort(mx.cbegin(), mx.cend(),
      [] (const auto & p1, const auto & p2)
      {
        return (p1.second > p2.second);
      }
    );

    if (mx[0].second == 1)
    {
      return mx;
    }
    else if (mx.size() < 5)
    {
      return mx;
    }

    int trun_num = std::max(2, int(floor(0.05 * mx.size())));
    float rate = 1/(
      // the average of the second values in each pair, cutting off at trunc_num
      std::accumulate(mx.begin(), mx.end() - trun_num, 0, 
        [] (const auto & p1, const auto & p2)
        {
          return (p1.second + p2.second);
        }
      ) / mx.size()
    );

    std::vector<int> 
    prob;

    for (const auto & i : mx)
    {
      prob.emplace_back(rate * exp(-rate * i.second));
    }
    
    // create and populate a cumulative_probability variable
    std::vector<int>
    cumulative_probability;
    int running_total;
    for (const auto & it : prob)
    {
      running_total += it;
      cumulative_probability.push_back(running_total);
    }

    if (cumulative_probability[0] > fdr_cutoff)
    {
      return mx;
    }
    else
    {
      std::vector<std::pair<int, int>>
      sliced_mx;

      for (const auto & i : cumulative_probability)
      {
        if (i > fdr_cutoff)
        {
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

    if (fs.size() > 1)
    {
      for (int it = 1; it != fs.size(); it++)
      {
        int
        clo_p = take_closest(tmp, fs[it]);

        if (abs(clo_p - fs[it]) > this->MAX_TS_DIST / 2)
        {
          tmp.push_back(fs[it]);
        }
      }
    }

    // if known_site is not None
    if (known_site.size() > 0)
    {
      for (auto s : known_site)
      {
        int
        clo_p = take_closest(tmp, s);
        if (abs(clo_p - s) > this->MAX_TS_DIST)
        {
          tmp.push_back(s);
        }
      }
    }

    return tmp;
  };
  
  // make counts objects for both left and right
  std::map<int, int>
  left_counts;

  for (const auto & i : this->left)
  {
    left_counts[i]++;
  }

  std::vector<int>
  cnt_l;

  std::map<int, int>
  right_counts;

  for (const auto & i : this->right)
  {
    right_counts[i]++;
  }

  std::vector<int>
  cnt_r;

  if ((left_counts.size() < this->MIN_SUP_CNT) || (right_counts.size() < this->MIN_SUP_CNT))
  {
    return;
  }
  else
  {
    // left
    std::vector<std::pair<int, int>>
    fs_l = filter_site(left_counts);

    
    if (fs_l.size() == 0)
    {
      cnt_l = {-99999999};
    }
    else
    {
      // print everything to out_f
      for (const auto & it : fs_l)
      {
        out_f << this->ch << "\t" 
              << it.first << "\t" 
              << it.first + 1 << "\t" 
              << it.second << "\n";
      }

      // make a vector to store the keys from fs_l
      std::vector<int> 
      fs_l_keys;
      for (const auto & i : fs_l)
      {
        fs_l_keys.push_back(i.first);
      }

      cnt_l = insert_dist(fs_l_keys, known_site.left);
      std::sort(cnt_l.cbegin(), cnt_l.cend());
    }

    // right
    std::vector<std::pair<int, int>>
    fs_r = filter_site(right_counts);

    if (fs_r.size() == 0)
    {
      cnt_r = {-99999999};
    }
    else
    {
      // print everything to out_f
      for (const auto & it : fs_r)
      {
        out_f << this->ch << "\t"
              << it.first << "\t"
              << it.first + 1 << "\t"
              << it.second << "\n";
      }

      // make a vector to store the keys from fs_r
      std::vector<int>
      fs_r_keys;
      for (const auto & i : fs_r)
      {
        fs_r_keys.push_back(i.first);
      }

      cnt_r = insert_dist(fs_r_keys, known_site.right);
      std::sort(cnt_r.cbegin(), cnt_r.cend());
    }
  }

  for (const auto & pair : lr_pair)
  {
    // first we need a counts map, which we populate
    std::map<std::pair<int, int>, int>
    pair_counts;
    for (const auto & i : pair.second)
    {
      pair_counts[i] ++;
    }

    // then to sort it, we need to convert to a vector
    std::vector<std::pair<std::pair<int, int>, int>>
    tmp_pair;
    for (const auto & i : pair_counts)
    {
      tmp_pair.push_back(i);
    }
    std::sort(tmp_pair.cbegin(), tmp_pair.cend(),
      [] (const auto & p1, const auto & p2)
      {
        return (p1.second > p2.second);
      }
    );

    std::vector<std::pair<int, int>>
    pair_after_filtering;

    // lines like these make me think, 
    // maybe i've not thought this architecture out well enough...
    std::vector<std::pair<std::pair<int, int>, int>>
    pair_enrich;

    for (const auto & p : tmp_pair)
    {
      std::pair<int, int>
      cl_p = {
        take_closest(cnt_l, p.first.first),
        take_closest(cnt_r, p.first.second)
      };

      if ((abs(cl_p.first - p.first.first) < (0.5 * this->MAX_TS_DIST)) and
          (abs(cl_p.second - p.first.second) < 0.5 * this->MAX_TS_DIST))
      {
        if (pair_after_filtering.size() > 0)
        {
          // line 650 of sc_longread.py

          // now we want to extract the left and right elements
          std::vector<int>
          pair_after_filtering_left;
          std::vector<int>
          pair_after_filtering_right;

          for (const auto & it : pair_after_filtering)
          {
            pair_after_filtering_left.push_back(it.first);
            pair_after_filtering_right.push_back(it.second);
          }

          int distance = (
            abs(take_closest(pair_after_filtering_left, cl_p.first) - cl_p.first) +
            abs(take_closest(pair_after_filtering_right, cl_p.second) - cl_p.second)
          );

          if (distance > this->MAX_TS_DIST)
          {
            if ((cl_p.first < pair.first.front()) and cl_p.second > pair.first.back())
            {
              pair_after_filtering.push_back(cl_p);
            }
          }
        }
        else
        {
          if ((cl_p.first < pair.first.front()) and (cl_p.second > pair.first.back()))
          {
            pair_after_filtering.push_back(cl_p);
          }
        }

        // but we only want to allow up to MAX_SITE_PER_SLICE combinations
        if (pair_after_filtering.size() >= this->MAX_SITE_PER_SLICE)
        {
          break;
        }
      }
    }

    if (pair_after_filtering.size() == 0)
    {
      // then search for isoform-specific enrichment
      for (const auto & p : tmp_pair)
      {
        // we need to sum up all of the elements in tmp_pair that meet our criteria
        int sum = 0;
        for (const auto & it : tmp_pair)
        {
          if ((abs(it.first.first - p.first.first) + abs(it.first.second - p.first.second)) < this->MAX_TS_DIST)
          {
            sum += it.second;
          }
        }
        pair_enrich.push_back({p.first, sum});
      }
      std::sort(
        pair_enrich.cbegin(), pair_enrich.cend(),
        [] (const auto & p1, const auto & p2)
        {
          // sign is flipped because we are reverse sorting
          return (p1.second < p2.second);
        }
      );

      for (const auto & p : pair_enrich)
      {
        if (pair_after_filtering.size() > 0)
        {
          // now we want to extract the left and right elements again
          std::vector<int>
          pair_after_filtering_left;
          std::vector<int>
          pair_after_filtering_right;

          for (const auto & it : pair_after_filtering)
          {
            pair_after_filtering_left.push_back(it.first);
            pair_after_filtering_right.push_back(it.second);
          }

          int distance = (
            abs(take_closest(pair_after_filtering_left, p.first.first) - p.first.first) +
            abs(take_closest(pair_after_filtering_right, p.first.second) - p.first.second)
          );

          if ((distance > this->MAX_TS_DIST) and 
              (p.first.first < pair.first.front()) and 
              (p.first.second > pair.first.back()))
          {
            // append the key of p
            pair_after_filtering.push_back(p.first);
          }
        }
        else if ((p.first.first < pair.first.front()) and (p.first.second > pair.first.back()))
        {
          // append the key of p
          pair_after_filtering.push_back(p.first);
        }

        // we only want up to MAX_SITE_PER_SLICE combinations
        if (pair_after_filtering.size() >= this->MAX_SITE_PER_SLICE)
        {
          // so just disreagard any after that point
          break;
        }
      }
    }

    int sup_cnt_total = 0;

    for (const auto & p : pair_after_filtering)
    {
      // now we want to add the values of keys matching certain criteria
      for (const auto & it : tmp_pair)
      {
        if ((abs(it.first.first - p.first) + abs(it.first.second - p.second)) < this->MAX_TS_DIST)
        {
          sup_cnt_total += it.second;
        }
      }
    }

    if ((float(sup_cnt_total) / this->lr_pair.size()) <= this->MIN_SUP_PCT)
    {
      // not enough support counts - this often happens in the case when the TSS/TES is strongly degraded

      if (pair_enrich.size() == 0)
      {
        for (const auto & p : tmp_pair)
        {
          // add it to pair_enrich
          // we need to sum up all of the elements in tmp_pair that meet our criteria
          int sum = 0;
          for (const auto & it : tmp_pair)
          {
            if ((abs(it.first.first - p.first.first) + abs(it.first.second - p.first.second)) < this->MAX_TS_DIST)
            {
              sum += it.second;
            }
          }
          pair_enrich.push_back({p.first, sum});
        }
        std::sort(
          pair_enrich.cbegin(), pair_enrich.cend(),
          [] (const auto & p1, const auto & p2)
          {
            // sign is flipped because we are reverse sorting
            return (p1.second < p2.second);
          }
        );
      }

      std::vector<int>
      tmp_ex;
      tmp_ex.push_back(int(pair_enrich[0].first.first));
      for (int i : pair.first)
      {
        tmp_ex.push_back(i);
      }
      tmp_ex.push_back(int(pair_enrich[0].first.second));

      this->raw_isoforms[tmp_ex] = pair.second;
      this->strand_cnt[tmp_ex] = this->strand_cnt[pair.first];
    }
    else
    {
      // now, add filtered TSS/TES to raw_isoforms
      for (const auto & p : pair_after_filtering)
      {
        std::vector<int>
        tmp_ex;
        tmp_ex.push_back(int(p.first));
        for (int i : pair.first)
        {
          tmp_ex.push_back(i);
        }
        tmp_ex.push_back(int(p.second));

        // we need to sum up all of the elements in tmp_pair that meet our criteria
        int sum = 0;
        for (const auto & it : tmp_pair)
        {
          if ((abs(it.first.first - p.first) + abs(it.first.second - p.second)) < this->MAX_TS_DIST)
          {
            sum += it.second;
          }
        }

        // come back and change this line later
        this->raw_isoforms[tmp_ex] = {{sum, 0}};
        this->strand_cnt[tmp_ex] = this->strand_cnt[pair.first];
      }
    }
  }
}

void Isoforms::match_known_annotation (
  std::map<std::string, Junctions> transcript_to_junctions,
  std::map<std::string, Pos> transcript_dict,
  std::map<std::string, int> gene_dict,
  GeneBlocks one_block,
  std::map<std::string, int> fa_dict
)
{

  std::set<int>
  splice_site = get_splice_site(transcript_to_junctions, one_block.transcript_list);

  std::vector<std::vector<int>>
  junction_list;
  std::map<std::vector<int>, std::string>
  junction_dict;

  std::vector<std::vector<int>>
  exons_list;
  std::map<std::vector<int>, std::string>
  exons_dict;

  // populate junction_list, junction_dictionary
  for (const auto & i : one_block.transcript_list)
  {
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

  if (splice_site.size() == 0)
  {
    // if this is a single-exon read, we're done
    return;
  }

  for (const auto & [exon_key, exon_val] : this->single_block_dict)
  {
    for (const auto & i : one_block.transcript_list)
    {
      std::vector<int>
      tmp_std;

      if (this->STRAND_SPECIFIC == 0)
      {
        // we need to convert the pair to a vector for this lookup
        tmp_std = this->strand_cnt[{exon_key.first, exon_key.second}];
      }
      else
      {
        if (transcript_dict[i].strand == '+')
        {
          tmp_std = {1};
        }
        else
        {
          tmp_std = {-1};
        }
      }

      if ((transcript_to_junctions[i].junctions.size() == 0) &&
          (tmp_std == this->strand_cnt[{exon_key.first, exon_key.second}]))
      {
        if ((abs(exon_key.first - transcript_to_junctions[i].left[0]) < this->MAX_TS_DIST) &&
            (abs(exon_key.second - transcript_to_junctions[i].right[0]) < this->MAX_TS_DIST))
        {
          std::pair<int, int>
          known_exons = {transcript_to_junctions[i].left[0], transcript_to_junctions[i].right[0]};
          if (this->known_isoforms.count({known_exons.first, known_exons.second})) // check if it is already known
          {
            // it's already known, just update the entry
            this->known_isoforms[{known_exons.first, known_exons.second}] = Iso(
              this->known_isoforms[{known_exons.first, known_exons.second}].support_cnt + this->single_blocks[this->single_block_dict[exon_key]].size(),
              i,
              transcript_dict[i].parent_id
            );
          }
          else
          {
            // it's totally new, create a fresh entry for it
            this->known_isoforms[{known_exons.first, known_exons.second}] = Iso(
              this->single_blocks[this->single_block_dict[exon_key]].size(),
              i,
              transcript_dict[i].parent_id
            );
            if (!this->STRAND_SPECIFIC)
            {
              this->strand_cnt[{known_exons.first, known_exons.second}] = {1};
              if (transcript_dict[i].strand == '-')
              {
                this->strand_cnt[{known_exons.first, known_exons.second}] = {-1};
              }
            }
          }
        }
      }
    }
  }

  for (auto const & [raw_iso_key, raw_iso_val] : this->raw_isoforms)
  {
    // line 730 in sc_longread.py
    int found = 0;

    // we need to check if there are any repeated values in raw_iso_key
    std::set<int>
    raw_iso_key_set;
    for (const auto & i : raw_iso_key)
    {
      raw_iso_key_set.insert(i);
    }
    if (raw_iso_key.size() > raw_iso_key_set.size())
    {
      // ignore raw_iso data that has repeated entries
      continue;
    }

    std::vector<int>
    tmp_std;

    for (const auto & i : one_block.transcript_list)
    {
      if (!this->STRAND_SPECIFIC)
      {
        tmp_std = this->strand_cnt[raw_iso_key];
      }
      else
      {
        tmp_std = {1};
        if (transcript_dict[i].strand == '+')
        {
          tmp_std = {-1};
        }
      }

      // same number of exons, same strand
      if ((raw_iso_key.size() - 2 == transcript_to_junctions[i].junctions.size()) &&
          (tmp_std == this->strand_cnt[raw_iso_key]))
      {
        int iso_is_within_max_dist = 1;
        for (int j = 0; j < std::min(raw_iso_key.size() - 1, transcript_to_junctions[i].junctions.size()); j++)
        {
          if (abs(raw_iso_key[j+1] - transcript_to_junctions[i].junctions[j]) > this->MAX_DIST)
          {
            // we don't want any to be greater than MAX_DIST
            iso_is_within_max_dist = 0;
            break;
          }
        }

        if (iso_is_within_max_dist)
        {
          if ((abs(raw_iso_key.front() - transcript_to_junctions[i].left[0]) < this->MAX_DIST) &&
              (abs(raw_iso_key.back() - transcript_to_junctions[i].right[0]) < this->MAX_DIST))
          {
            // populate known_exons with transcript_to_junctions
            std::vector<int>
            known_exons = {transcript_to_junctions[i].left[0]};
            for (const auto & j : transcript_to_junctions[i].junctions)
            {
              known_exons.push_back(j);
            }
            known_exons.push_back(transcript_to_junctions[i].right[0]);

            if (this->known_isoforms.count(known_exons))
            {
              // line 747 sc_longread.py
              
          }
        }
      }
    }
  }
}