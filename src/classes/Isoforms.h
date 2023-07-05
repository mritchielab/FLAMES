#ifndef ISOFORMS
#define ISOFORMS

#include <map>
#include <unordered_map>
#include <vector>
#include <tuple>

#include <Rcpp.h>

#include "StartEndPair.h"
#include "Pos.h"
#include "junctions.h"
#include "GeneBlocks.h"

/*
    a data container used in Isoforms,
    specifically for known_isoforms and match_known_annotation
*/
struct Iso
{
	long support_count;
	std::string transcript_id; 
	std::string gene_id;

	bool operator==(const Iso &rhs) {
		return std::tie(support_count, transcript_id, gene_id) == std::tie(rhs.support_count, rhs.transcript_id, rhs.gene_id);
	}
	bool operator!=(const Iso &rhs) {
		return std::tie(support_count, transcript_id, gene_id) != std::tie(rhs.support_count, rhs.transcript_id, rhs.gene_id);
	}
};

class Isoforms {
private:
    // these values will all be extracted from config
    // Rcpp::List parameters;
    int MAX_TS_DIST;
    int MAX_DIST;
    int strand_specific;
    int MAX_SPLICE_MATCH_DIST;
    int remove_incomp_reads;
    int min_fl_exon_len;
    int Min_sup_cnt;
    float Min_sup_pct;
    int Max_site_per_splice;

    void add_one(Junctions junctions, bool strand);
    void update_one(Junctions junctions, std::vector<int> key, bool strand);
    void update_one(Junctions junctions, StartEndPair key, bool strand);

    void update_new_isoform(const std::vector<int> &, long, const std::string &, const std::string &);

    std::vector<int> insert_dist(std::vector<int>, std::vector<int>);

    static std::vector<std::pair<int, int>> filter_site(const std::unordered_map<int, int> &, float);
    static inline void output_to_bedgraph(std::ofstream &, std::string, int, int, int);
public:
    std::string ch;
    // map of junctions to vectors of junctions
    std::unordered_map<std::vector<int>, std::vector<std::vector<int>>>
    junction_dict;
    std::unordered_map<std::vector<int>, std::vector<StartEndPair>>
    lr_pair;

    std::vector<int> 
    left;
    std::vector<int> 
    right;

    std::unordered_map<StartEndPair, std::vector<StartEndPair>>
    single_block_dict;

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

    Isoforms(const std::string _ch, int _MAX_TS_DIST, int _MAX_DIST, int _strand_specific, int _MAX_SPLICE_MATCH_DIST,
            int _remove_incomp_reads, int _min_fl_exon_len, int _Min_sup_cnt,
            float _Min_sup_pct, int _Max_site_per_splice)
        : MAX_TS_DIST{_MAX_TS_DIST}, MAX_DIST{_MAX_DIST}, strand_specific{_strand_specific}, 
        MAX_SPLICE_MATCH_DIST{_MAX_SPLICE_MATCH_DIST}, remove_incomp_reads{_remove_incomp_reads}, 
        min_fl_exon_len{_min_fl_exon_len}, Min_sup_cnt{_Min_sup_cnt}, Min_sup_pct{_Min_sup_pct}, 
        Max_site_per_splice{_Max_site_per_splice}, ch{_ch} {}
    
    int size();
    void update_all_splice();

    void add_isoform(Junctions junctions, bool is_reversed);

    void filter_TSS_TES(std::ofstream &out_f, DoubleJunctions known_site={}, float fdr_cutoff=0.01); 

    void match_known_annotation (
        const std::unordered_map<std::string, Junctions> 					& transcript_to_junctions,
        const std::unordered_map<std::string, Pos> 							& transcript_dict,
        const std::unordered_map<std::string, std::vector<StartEndPair>> 	& gene_dict,
        GeneBlocks one_block,
        std::unordered_map<std::string, std::string> fa_dict
    );

    std::string isoform_to_gff3(float isoform_pct);

    std::vector<int> countLR(std::ofstream &, const std::unordered_map<int, int> &, std::vector<int>, float);
};

#endif // ISOFORMS