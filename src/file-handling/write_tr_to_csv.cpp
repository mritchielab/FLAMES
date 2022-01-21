#include <fstream>
#include <unordered_set>
#include <unordered_map>
#include <map>
#include <vector>
#include <string>
#include <numeric>
#include <algorithm>
#include <Rcpp.h>

#include "../classes/Pos.h"
#include "write_tr_to_csv.h"


int
edit_distance(std::string s1, std::string s2)
{
    /*  
        calculates Levenshtein distance between two strings
        (code adapted from:
            https://www.programcreek.com/cpp/?code=salman-bhai%2FDS-Algo-Handbook%2FDS-Algo-Handbook-master%2FAlgorithms%2FDynamic_Programming%2FEdit_Distance%2FeditDistanceTopDown.cpp
        )
    */
    int m = s1.size();
    int n = s2.size();
    int dp[m+1][n+1];

	// Rcpp::Rcout << "Starting edit_distance\n";
    for (int i = 0; i <= m; i++) {
		// Rcpp::Rcout << "\ti=" << i << "\n";
        for (int j = 0; j <= n; j++) {
			// Rcpp::Rcout << "\t\tj=" << j << "\n";
            if (i == 0) {
                // if first string is empty,
                // minimum j operations in that row
                dp[i][j] = j;
            } else if (j == 0) {
                // if second string is empty,
                // minimum i operations for that row
                dp[i][j] = i;
            } else if (s1[i-1] == s2[j-1]) {
                // if last characters are the same, go ahead for the rest of the characters
                dp[i][j] = dp[i-1][j-1];
            } else {
                dp[i][j] = 1 + std::min(std::min(dp[i][j-1], dp[i-1][j]), dp[i-1][j-1]);
            }
        }
    }

    return dp[m][n];
}

// [[Rcpp::export]]
void xx() {
	std::vector<std::pair<std::string, std::string>> ps {
			{"AATGC", "AAGGC"},
			{"AAGGGGC", "AATGA"},
			{"GTCGAGATCGA", "AGATCGTAGC"},
			{"GGGGGGGTGGGCG", "GGGGGGGGGGGAAG"},
			{"AAAAA", "AAAAA"}
		};
	std::vector<int> res = {
		1,
		4,
		7,
		3,
		0
	};

	std::string x1, x2;
	for (int i = 0; i < res.size(); i++) {
		x1 = ps[i].first;
		x2 = ps[i].second;
		Rcpp::Rcout << "String 1 : " << x1 << " String 2 : " << x2 << "\n";

		Rcpp::Rcout << edit_distance(x1, x2) << "\n";
	}

}

int
umi_dedup (std::vector<std::string> l, bool has_UMI)
{

	if (!has_UMI) {
		return l.size();
	}

	std::map<std::string, int> l_count;
	for (const auto & i : l) {
		l_count[i]++;
	}

	if (l_count.size() == 1) {
		return 1;
	}

	// to sort it, we need to convert to a vector
	std::vector<std::pair<std::string, int>> l_count_vec;
	for (const auto & pair : l_count) {
		l_count_vec.push_back(pair);
	}
	std::sort(l_count_vec.begin(), l_count_vec.end(),
		[] (const std::pair<std::string, int> & p1, const std::pair<std::string, int> & p2) {
			return (p1.second > p2.second);
		}
	);

	
	std::unordered_map<std::string, int> rm_umi;
	for (int i = 0; i < l_count.size(); ++i) {
		for (int j = l_count.size(); j > i; j--) { // first assess the low abundant UMI
			if (!rm_umi.count(l_count_vec[j].first)) {
				if (edit_distance(l_count_vec[i].first, l_count_vec[j].first) < 2) {
					rm_umi[l_count_vec[j].first] = 1;
				}
			}
		}
	}

	return l_count.size() - rm_umi.size();
}

std::unordered_map<std::string, int>
write_tr_to_csv_cpp
(
    std::unordered_map<std::string, std::unordered_map<std::string, std::vector<std::string>>>
    bc_tr_count_dict,
    std::unordered_map<std::string, Pos>
    transcript_dict,
    std::string
    csv_f,
    std::unordered_map<std::string, Pos>
    transcript_dict_ref,
    bool
    has_UMI
)
{
    Rcpp::Rcout << "started write_tr_to_csv_cpp\n";

    std::ofstream csv (csv_f);

    std::unordered_set<std::string> all_tr = {};

    for (const auto & [bc, tr_dict] : bc_tr_count_dict) {
        for (const auto & [tr, entries] : tr_dict) {
			all_tr.insert(tr);
        }
    }

    Rcpp::Rcout << "all_tr is " << all_tr.size() << " long\n";

    // write the header to the csv
    csv << "transcript_id,gene_id";
    for (const auto & [bc, tr_dict] : bc_tr_count_dict) {
        csv << "," << bc;
    }
    csv << "\n";

    std::unordered_map<std::string, int> tr_count;
    // get a sum of the number of unique (non-duplicate) UMIs for each tr 
    for (const auto & tr : all_tr) {
        std::vector<int> count_l = {};
        for (auto & [bc, tr_count_dict] : bc_tr_count_dict) {
            if (tr_count_dict.count(tr) > 0) {
                count_l.push_back(umi_dedup(tr_count_dict[tr], has_UMI)); 
                tr_count[tr] += count_l.back();
            } else {
                count_l.push_back(0);
            }
        }

        if (transcript_dict.count(tr)) {
            csv << tr << "," << transcript_dict[tr].parent_id << ",";
        } else if (transcript_dict_ref.size() && transcript_dict_ref.count(tr)) {
            csv << tr << "," << transcript_dict_ref[tr].parent_id << ",";
        } else {
            Rcpp::Rcout << "Cannot find transcript in transcript_dict: " << tr << "\n";
            return {};
        }

		// write count_l to the file as a comma separated list
		csv << std::accumulate(count_l.begin()+1, count_l.end(), std::to_string(*count_l.begin()), 
			[](std::string &acc, int cur){
				return acc + "," + std::to_string(cur);
			});
		csv << "\n";
    }

    Rcpp::Rcout << "tr_count is " << tr_count.size() << " long\n";
    
    csv << "\n";
	csv.close();
    return tr_count;
}