#include <fstream>
#include <unordered_set>
#include <unordered_map>
#include <map>
#include <vector>
#include <string>

#include <Rcpp.h>

#include "../classes/Pos.h"

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
    transcript_dict_ref={},

    bool
    has_UMI=true
);