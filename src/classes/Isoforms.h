#ifndef ISOFORMS
#define ISOFORMS

#include <map>
#include <unordered_map>
#include <vector>

#include "StartEndPair.h"
#include "Pos.h"
#include "../utility/junctions.h"
#include "GeneBlocks.h"
#include "Config.h"

struct Iso
{
	/*
		a data container used in Isoforms,
		specifically for known_isoforms and match_known_annotation
	*/

	long support_count;
	std::string transcript_id; 
	std::string gene_id;
};

class Isoforms {
	private:
		// these values will all be extracted from config
		IsoformParameters parameters;

		void add_one(Junctions junctions, bool strand);
		void update_one(Junctions junctions, std::vector<int> key, bool strand);
		void update_one(Junctions junctions, StartEndPair key, bool strand);

		std::vector<int> insert_dist(std::vector<int>, std::vector<int>);

		static std::vector<std::pair<int, int>> filter_site(const std::unordered_map<int, int> &, float);

		static void output_to_bedgraph(std::ofstream &, std::string, int, int, int);

		std::vector<int> countLR(std::ofstream &, const std::unordered_map<int, int> &, std::vector<int>, float);
  	public:
		Isoforms(const std::string ch, IsoformParameters parameters)
		{
		/*
			initialises the object,
		*/
		this->parameters = parameters;
		this->ch = ch;
		}
		
		std::string ch;
		
		// std::unordered_map<std::vector<int>, int> 
		// junction_dict;
		// std::vector<std::vector<std::vector<int>>> 
		// junction_list;
		// map of junctions to vectors of junctions
		std::unordered_map<std::vector<int>, std::vector<std::vector<int>>>
		junction_dict;
		std::unordered_map<std::vector<int>, std::vector<StartEndPair>>
		lr_pair;
		std::vector<int> 
		left;
		std::vector<int> 
		right;

		// std::unordered_map<StartEndPair, int> 
		// single_block_dict = {};
		// std::vector<std::vector<StartEndPair>> 
		// single_blocks;
		std::unordered_map<StartEndPair, std::vector<StartEndPair>>
		single_block_dict;
		// std::unordered_map<std::vector<int>, int> 
		// strand_counts;
		std::unordered_map<std::vector<int>, std::vector<int>>
		strand_counts;
		std::unordered_map<std::vector<int>, Iso> 
		new_isoforms;
		std::unordered_map<std::vector<int>, Iso> 
		known_isoforms;
		std::unordered_map<std::vector<int>, int> 
		raw_isoforms;
		std::unordered_map<std::string, std::vector<std::vector<int>>> 
		ge_dict;

		void add_isoform(Junctions junctions, bool is_reversed);
		int size();
		void update_all_splice();

		void filter_TSS_TES(std::ofstream &out_f, DoubleJunctions known_site={}, float fdr_cutoff=0.01); 

		//unused
		// std::pair<std::vector<int>, std::unordered_map<int, int>>
		// group_sites(std::vector<int> l, int smooth_window, int min_threshold);

		void match_known_annotation (
			const std::unordered_map<std::string, Junctions> 					& transcript_to_junctions,
			const std::unordered_map<std::string, Pos> 						& transcript_dict,
			const std::unordered_map<std::string, std::vector<StartEndPair>> 	& gene_dict,
			GeneBlocks one_block,
			std::unordered_map<std::string, std::string> fa_dict
		);

		std::string isoform_to_gff3(float isoform_pct);
};

#endif // ISOFORMS