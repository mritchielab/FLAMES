#ifndef WRITE_TR_TO_CSV_H
#define WRITE_TR_TO_CSV_H

#include <unordered_map>
#include <vector>
#include <string>

#include "Pos.h"

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

int
umi_dedup (std::vector<std::string>, bool);

int
edit_distance(std::string, std::string);

#endif // WRITE_TR_TO_CSV_H