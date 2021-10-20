#include "parse_realigned_bam.hpp"

std::unordered_map<std::string, int>
file_to_map(std::string filename)
{
    /* 
        parse a file into a map. the file should look like:
            text 10
            word 4
            dhgukrahgk 10
    */
    std::unordered_map<std::string, int>
    map;

    std::ifstream
    file(filename);

    std::string line;
    while (std::getline(file, line)) {
        std::stringstream
        linestream (line);

        // break each line into words
        std::vector<std::string> words;
        std::string word;
        while (std::getline(linestream, word, ' ')) {
            words.push_back(word);
        }

        map[words[0]] = atoi(words[1].c_str());
    }

    return map;
}

void
parse_realigned_bam(
    std::string bam_in,
    std::string fa_idx_f,
    std::string min_sup_reads,
    std::string min_tr_coverage,
    std::string min_read_coverage,
    std::string kwargs
)
{
    // we need to read in the fa_idx_f file line by line, adding each one to the dict
    std::unordered_map<std::string, int>
    fa_idx = file_to_map(fa_idx_f);

    std::unordered_map<std::string, std::string>
    bc_tr_count_dict = {};

    std::unordered_map<std::string, std::string>
    bc_tr_badvoc_count_dict = {};

    std::unordered_map<std::string, std::string>
    tr_cov_dict = {};
    
    std::unordered_map<std::string, int>
    count_stat = {};

    // bamfile
}