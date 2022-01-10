#ifndef ISOFORMS
#define ISOFORMS

#include <map>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <sstream>

#include <iostream>
#include <fstream>

#include "StartEndPair.h"
#include "Pos.h"
#include "../junctions.h"
#include "../utilities/misc.h"
#include "GeneBlocks.h"
#include "Config.h"

class Isoforms
{
  private:
    // these values will all be extracted from config
    IsoformParameters parameters;
    
  public:
    Isoforms(std::string ch, IsoformParameters parameters)
    {
      /*
        initialises the object,
        extracting all the const values from config

        (why are config keys a mix of CAPS, lower, and Sentence Case?)
        (I have no idea, that's how they were when I found them)
        (change it later)
      */
      this->parameters = parameters;
      this->ch = ch;
    }
    
    std::string ch;
    
    std::unordered_map<std::vector<int>, int> 
    junction_dict = {};
    std::vector<std::vector<std::vector<int>>> 
    junction_list;

    std::unordered_map<std::vector<int>, std::vector<StartEndPair>>
    lr_pair;
    std::vector<int> 
    left;
    std::vector<int> 
    right;

    std::unordered_map<StartEndPair, int> 
    single_block_dict = {};
    std::vector<std::vector<StartEndPair>> 
    single_blocks;
    std::unordered_map<std::vector<int>, int> 
    strand_counts;
    std::unordered_map<std::vector<int>, std::vector<int>>
    new_strand_counts;
    std::unordered_map<std::vector<int>, Iso> 
    new_isoforms;
    std::unordered_map<std::vector<int>, Iso> 
    known_isoforms;
    std::unordered_map<std::vector<int>, int> 
    raw_isoforms;
    std::unordered_map<std::string, std::vector<std::vector<int>>> 
    ge_dict;

    void add_isoform(Junctions junctions, bool is_reversed);
    void add_one(Junctions junctions, bool strand);
    void update_one(Junctions junctions, std::vector<int> key, bool strand);
    int len();
    void update_all_splice();

    void filter_TSS_TES(std::ofstream * out_f, Junctions known_site={}, float fdr_cutoff=0.01); 

    //unused
    std::pair<std::vector<int>, std::unordered_map<int, int>>
    group_sites(std::vector<int> l, int smooth_window, int min_threshold);

    void match_known_annotation (
      std::unordered_map<std::string, Junctions> transcript_to_junctions,
      std::unordered_map<std::string, Pos> transcript_dict,
      std::unordered_map<std::string, std::vector<StartEndPair>> gene_dict,
      GeneBlocks one_block,
      std::unordered_map<std::string, std::string> fa_dict
    );

    std::string isoform_to_gff3(float isoform_pct);
};

#endif