#ifndef GFF3_TO_FA_H
#define GFF3_TO_FA_H

#include <string>
#include <vector>
#include <unordered_map>
#include <fstream>
#include <sstream>

#include "../classes/Pos.h"
#include "../classes/ReferenceDict.h"
#include "../classes/StartEndPair.h"

void
get_transcript_seq
(
    std::string fa_file,
    std::string fa_out_f,
    std::unordered_map<std::string, std::vector<std::string>>   * chr_to_gene,
    std::unordered_map<std::string, Pos>                        * transcript_dict,
    std::unordered_map<std::string, std::vector<std::string>>   * gene_to_transcript,
    std::unordered_map<std::string, std::vector<StartEndPair>>  * transcript_to_exon,

    ReferenceDict * ref_dict = nullptr
);


std::unordered_map<std::string, std::string>
get_fa_simple(std::string filename);

void
write_fa(std::ofstream* fa_out, std::string na, std::string seq, int wrap_len=50);

std::string
r_c(const std::string * seq);

#endif // GFF3_TO_FA_H