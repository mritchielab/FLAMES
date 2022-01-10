#include "write_tr_to_csv.hpp"


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
    int dp[m][n];

    for (int i = 0; i <= m; ++i) {
        for (int j = 0; j <= n; j++) {
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

int
umi_dedup
(
    std::vector<std::string> l, 
    bool has_UMI
)
{
    if (has_UMI) {
        std::map<std::string, int>
        l_count;
        for (const auto & i : l) {
            if (l_count.count(i) == 0) {
                l_count[i] = 0;
            }
            l_count[i]++;
        }

        if (l_count.size() == 1) {
            return 1;
        }

        // to sort it, we need to convert to a vector
        std::vector<std::pair<std::string, int>>
        l_count_vec;
        for (const auto & pair : l_count) {
            l_count_vec.push_back(pair);
        }
        std::sort(l_count_vec.begin(), l_count_vec.end(),
            [] (const auto & p1, const auto & p2) {
                return (p1.second > p2.second);
            }
        );

        
        std::vector<std::string>
        rm_umi;

        for (int i = 0; i < l_count.size(); ++i) {
            for (int j = l_count.size(); j > i; j--) { // first assess the low abundant UMI
                if (std::count(rm_umi.begin(), rm_umi.end(), l_count_vec[j].first) == 0) { 
                    if (edit_distance(l_count_vec[i].first, l_count_vec[j].first) < 2) {
                        rm_umi.push_back(l_count_vec[j].first);
                    }
                }
            }
        }

        return l_count.size() - rm_umi.size();
    } else {
        return l.size();
    }
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
    std::cout << "started write_tr_to_csv_cpp\n";

    std::cout << "bc_tr_count_dict: (size "<< bc_tr_count_dict.size() << ")\n";
    for (const auto & [bc, tr] : bc_tr_count_dict) {
        std::cout << bc << ": (size " << tr.size() << ") ";
        for (const auto & st : tr) {
            std::cout << st.first << ", ";
        }
    }

    std::cout << "transcript_dict: (size " << transcript_dict.size() << ")\n";
    for (const auto & [tr, pos] : transcript_dict) {
        std::cout << tr << ":(" << pos.start << "," << pos.end << ")\n";
    }

    std::cout << "transcript_dict_ref: (size " << transcript_dict_ref.size() << ")\n";
    for (const auto & [tr, pos] : transcript_dict_ref) {
        std::cout << tr << ":("<<pos.start<<","<<pos.end<<")\n";
    }

    std::ofstream
    csv (csv_f);

    std::unordered_set<std::string> 
    all_tr = {};

    for (const auto & [bc, tr_dict] : bc_tr_count_dict) {
        for (const auto & [tr, entries] : tr_dict) {
            if (all_tr.count(tr) == 0) {
                all_tr.insert(tr);
            }
        }
    }

    std::cout << "all_tr is " << all_tr.size() << " long\n";

    // write the header to the csv
    csv << "transcript_id,gene_id";
    for (const auto & [bc, tr_dict] : bc_tr_count_dict) {
        csv << "," << bc;
    }
    csv << "\n";

    std::unordered_map<std::string, int>
    tr_count;
    // get a sum of the number of unique (non-duplicate) UMIs for each tr 
    for (const auto & tr : all_tr) {
        tr_count[tr] = 0;
        std::vector<int>
        count_l = {};
        for (auto & [bc, tr_count_dict] : bc_tr_count_dict) {
            if (tr_count_dict.count(tr) > 0) {
                count_l.push_back(umi_dedup(tr_count_dict[tr], has_UMI));
                tr_count[tr] += count_l.back();
            } else {
                count_l.push_back(0);
            }
        }

        if (transcript_dict.count(tr) > 0) {
            csv << tr << "," << transcript_dict[tr].parent_id << ",";
        } else if (transcript_dict_ref.size() > 0 && transcript_dict_ref.count(tr) > 0) {
            csv << tr << "," << transcript_dict_ref[tr].parent_id << ",";
        } else {
            std::cout << "Cannot find transcript in transcript_dict: " << tr << "\n";
            return {};
        }

        for (const auto & x : count_l) {
            csv << x << ",";
        }
    }

    std::cout << "tr_count is " << tr_count.size() << " long\n";
    
    csv << "\n";
    return tr_count;
}