#include "GFFData.h"

void
GFFData::remove_transcript_duplicates(bool update_transcript_dict) {
    // Remove duplicates from the transcript_to_exon maps
    for (auto tr : transcript_to_exon) {
        // it->first is std::string key
        // it->second is std::vector<StartEndPair> list of pairs
        std::sort(tr.second.begin(), tr.second.end(), StartEndPairCompare);
        // if the forward list has at least two elements and start of first element equals start of second
        if (tr.second.size() >= 2 &&
            (tr.second.begin()->start) == (tr.second.begin() + 1)->start) {
            
            std::vector<StartEndPair> new_ex = {*tr.second.begin()};
			// remove duplicates from the old list of exons stored in transcript_to_exon[x]
            for (auto ex : tr.second) {
				StartEndPair last_ex = *(new_ex.end() - 1);
				if (ex.start != last_ex.start && ex.end != last_ex.start) {
					new_ex.push_back(ex);
				}
            }

            transcript_to_exon[tr.first] = new_ex;
        }

		if (update_transcript_dict) {
                // update transcript_dict[this transcript] with a new Pos object of
                // the correct start and end positions of this exon
                transcript_dict[tr.first] = Pos {
                    transcript_dict[tr.first].chr,
                    tr.second.begin()->start,
                    (tr.second.end() - 1)->end,
                    transcript_dict[tr.first].strand,
                    transcript_dict[tr.first].parent_id
                };
            }
    }

	// remove duplicates from gene_to_transcript
	for (auto ge : gene_to_transcript) {
		std::unordered_map<std::string, int> set;
		std::vector<std::string> new_genes;

		for (auto el : ge.second) {
			set[el] = 1;
		}

		for (auto unique : set) {
			new_genes.push_back(unique.first);
		}

		ge.second = new_genes;
	}
}


List
GFFData::to_R()
{
    // /// Create a Rcpp::List from an unordered_map.
    // /// Specificially used for chr_to_gene and gene_to_transcript in order to export each object
    // /// to the R session for later use.
    // /// Maintains order of data
    // auto create_list_map_to_map = [] (std::unordered_map<std::string, std::unordered_map<std::string, bool>> map)
    // {
    //     List result = List::create();

    //     for (auto it = map.begin; it != map.end(); ++it) {
    //         // iterating over every string to unordered_map in map
    //         CharacterVector inner_list = CharacterVector::create();

    //         // iterating over every string in the internal unordered_map
    //         for (auto in = it->second.begin(); in != it->second.end(); ++in) {
    //             inner_list.push_back(in->first);
    //         }

    //         result.push_back(inner_list, it->first);
    //     }

    //     return result;
    // };

    // auto create_list_map_to_pos = [] (std::unordered_map<std::string, Pos> map)
    // {
    //     List result = List::create();

    //     // result should look like: "transcript" = ["pos"=[]]
    //     for (auto it = map.begin(); it != map.end(); ++it) {
    //         Pos pos = it->second;

    //         result.push_back(pos_to_R(&pos), it->first);
    //     }

    //     return result;
    // };

    // auto create_list_map_to_startendpair = [] (std::unordered_map<std::string, std::list<StartEndPair>> map)
    // {
    //     List transcript_to_exon_list = List::create();

    //     for (auto tr = map.begin(); tr !=map.end(); ++tr) {
    //         // it->first is std::string key
    //         // it->second is std::list<StartEndPair> list of pairs
    //         List pair_list = List::create();

    //         for (auto ex = tr->second.begin(); ex != tr->second.end(); ++ex) {
    //             pair_list.push_back(
    //                 IntegerVector::create(
    //                     ex->start,
    //                     ex->end
    //                 )
    //             );
    //         }

    //         transcript_to_exon_list.push_back(pair_list, tr->first);
    //     }

    //     return transcript_to_exon_list;
    // };

    // return List::create(
    //     _["chr_to_gene"] = create_list_map_to_map(this->chr_to_gene),
    //     _["transcript_dict"] = create_list_map_to_pos(this->transcript_dict),
    //     _["gene_to_transcript"] = create_list_map_to_map(this->gene_to_transcript),
    //     _["transcript_to_exon"] = create_list_map_to_startendpair(this->transcript_to_exon)
    // );
    return List::create();
}

void
GFFData::from_R(List list)
{
    /* import an R list */
}