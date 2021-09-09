#include <map>
#include <vector>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <sstream>

#include <iostream>
#include <fstream>

#include "junctions.h"
#include "misc.h"
#include "GeneBlocks.h"

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
    std::map<std::vector<int>, int> 
    strand_counts;
    std::map<std::vector<int>, std::vector<int>>
    new_strand_counts;
    std::map<std::vector<int>, Iso> 
    new_isoforms;
    std::map<std::vector<int>, Iso> 
    known_isoforms;
    std::map<std::vector<int>, int> 
    raw_isoforms;
    std::map<std::string, std::vector<std::vector<int>>> 
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

    void match_known_annotation (
      std::map<std::string, Junctions> transcript_to_junctions,
      std::map<std::string, Pos> transcript_dict,
      std::map<std::string, std::vector<int>> gene_dict,
      GeneBlocks one_block,
      std::map<std::string, std::vector<char>> fa_dict
    );

    std::string isoform_to_gtt3(int isoform_pct);
};