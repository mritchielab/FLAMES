#include "write_tr_to_csv.hpp"

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

        if (l_count.size() == 1) {
            return 1;
        }
        
        std::map<char, std::string>
        rm_umi;

        for (int i = 0; i < l_count.size(); ++i) {
            for (int j = l_count.size(); j > i; j--) {
                if (rm_umi.count(l_count[j][0]) == 0) {

                }
            }
        }
    }
}

void
write_tr_to_csv
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
    std::ofstream
    csv (csv_f);

    std::unordered_set<std::string> all_tr = {};

    for (const auto & [bc, tr_dict] : bc_tr_count_dict) {
        for (const auto & [tr, entries] : tr_dict) {
            if (all_tr.count(tr) == 0) {
                all_tr.insert(tr);
            }
        }
    }

    // write the header to the csv
    csv << "transcript_id,gene_id";
    for (const auto & [bc, tr_dict] : bc_tr_count_dict) {
        csv << "," << bc;
    }
    csv << "\n";

    std::unordered_map<std::string, int>
    tr_count;

    for (const auto & tr : all_tr) {
        int
        count_l = 0;
        for ()
    }
}