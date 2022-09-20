#include "parse_realigned_bam.h"

#include <unordered_map>
#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include <iostream>

#include "../classes/BamRecord.h"

#include "../utility/bam.h"
#include "../utility/utility.h"
#include "../main-functions/group_bam2isoform.h"

std::unordered_map<std::string, int>
file_to_map(std::string filename)
{
    /* 
        parse a file into a map. the file should look like:
            text 10
            word 4
            dhgukrahgk 10
    */
    std::unordered_map<std::string, int> map;

    std::ifstream file(filename);

    std::string line;
    int lines = 0;
    while (std::getline(file, line)) {
        lines++;
        
		std::pair<std::string, std::string> parsed = parseSpace(line);
		std::pair<std::string, std::string> intParsed = parseSpace(parsed.second);

		map[parsed.first] = stoi(intParsed.first);
    }

    return map;
}

// [[Rcpp::export]]
void
read_entire_bam
(
    std::string bam_in, std::string log_out
)
{
    // open the outfile
    std::ofstream
    log (log_out);

    // read a bamfile
    bamFile bam = bam_open(bam_in.c_str(), "r"); // bam.h
    // bam_index_t *bam_index = bam_index_load(bam_in.c_str());
    bam_header_t *header = bam_header_read(bam); // bam.h

    bam1_t *b = bam_init1();

    while (bam_read1(bam, b) >= 0) {
        // skip the invalid tid entries
        if (b->core.tid == -1) {
            continue;
        }
        BAMRecord rec = read_record(b, header);
        log << "tid=" << b->core.tid 
            << "\tflag=" << b->core.flag 
            << "\t\tname=" << rec.reference_name 
            << "\treference_start=" << rec.reference_start 
            << "\treference_end=" << rec.reference_end 
            << "\tread_name=" << rec.read_name
            << "\tquality=" << rec.mapping_quality << "\n";
        
        std::cout << rec.reference_name << "\n";
    }
    bam_close(bam);
}

RealignedBamData
parse_realigned_bam
(
    std::string bam_in,
    std::string fa_idx_f,
    int         min_sup_reads,
    float       min_tr_coverage,
    float       min_read_coverage,
    std::unordered_map<std::string, std::string> kwargs
)
{

    // we need to read in the fa_idx_f file line by line, adding each one to the dict
    std::unordered_map<std::string, int>
    fa_idx = file_to_map(fa_idx_f);

    std::unordered_map<std::string, std::unordered_map<std::string, std::vector<std::string>>>
    bc_tr_count_dict = {};

    std::unordered_map<std::string, std::unordered_map<std::string, std::vector<std::string>>>
    bc_tr_badcov_count_dict = {};

    std::unordered_map<std::string, std::vector<float>>
    tr_cov_dict = {};

    std::unordered_map<std::string, std::vector<ReadDictEntry>>
    read_dict;
    
    // std::unordered_map<std::string, int>
    // count_stat = {};
	// {
	// 	unmapped: 0, 
	// 	not enough coverage: 1, 
	// 	no good match: 2, 
	// 	counted reads: 3, 
	// 	ambiguous reads: 4, 
	// 	not in annotation: 5
	// }
	std::vector<int> count_stat = {0, 0, 0, 0, 0, 0};

    std::unordered_map<std::string, std::string>
    bc_dict = {};
    if (kwargs.count("bc_file") > 0) {
        bc_dict = make_bc_dict(kwargs["bc_file"]);
    }
    
    // read a bamfile
    bamFile bam = bam_open(bam_in.c_str(), "r"); // bam.h
    //bam_index_t *bam_index = bam_index_load(bam_in.c_str());
    bam_header_t *header = bam_header_read(bam); // bam.h

    // iterate over every entry in the bam
    bam1_t *b = bam_init1();
    while (bam_read1(bam, b) >= 0) {
        if (b->core.tid == -1) {
            continue;
        }
        BAMRecord rec = read_record(b, header);
    
        // if it's unmapped, just update the count and continue
        if (rec.flag.read_unmapped) {
            count_stat[0] += 1;
            continue;
        }

        int map_start = rec.reference_start;
        int map_end = rec.reference_end;

        std::string tr = rec.reference_name;

        float tr_cov = float(map_end - map_start) / float(fa_idx[tr]);

        // if there is no dictionary entry for this tr_cov, create one
        // then add the new cov
        tr_cov_dict[tr].push_back(tr_cov);

        float inferred_read_length = calculate_length_from_cigar(rec.cigar, true);
		float tr_length = float(rec.query_alignment_length) / inferred_read_length;

        // set up the entry we're about to add to the dictionary
        ReadDictEntry new_read_entry = {
            tr,
            rec.AS_tag,
            tr_cov,
            tr_length,
            rec.mapping_quality
        };

        // create a dictionary entry if there isn't already one
        if (read_dict.count(rec.read_name) == 0) {
            read_dict[rec.read_name].push_back(new_read_entry);
        } else {
            if (rec.AS_tag > read_dict[rec.read_name][0].AS_tag) {
                read_dict[rec.read_name].insert(read_dict[rec.read_name].begin(), new_read_entry);
            } else if (new_read_entry.AS_tag == read_dict[rec.read_name][0].AS_tag &&
                        new_read_entry.length == read_dict[rec.read_name][0].length) {
                if (new_read_entry.tr_cov > read_dict[rec.read_name][0].tr_cov) { // choose the one with higher transcript coverage, might be internal TSS
                    read_dict[rec.read_name].insert(read_dict[rec.read_name].begin(), new_read_entry);
                }
            } else {
                read_dict[rec.read_name].push_back(new_read_entry);
            }
        }

        // see if tr is in the fa_idx, if not then log it as "not in annotation"
        if (fa_idx.count(tr) == 0) {
            count_stat[5] += 1;
        }
    }

	bam_close(bam);

    // build up a new dict of only the tr entries that have sufficient coverage
    // (at least 90% coverage on at least min_sup_reads)
    std::vector<std::string> tr_kept;
    for (const auto & [tr, covs] : tr_cov_dict) {
		// count up how many of the coverages are above the acceptable value of 90%
        int sup_read_count = 0;
        for (const auto & cov : covs) {
            if (cov > 0.9f) {
                sup_read_count++;
            }
        }
        
        // if the overall threshold is met for this tr, add it to the new dict
        if (sup_read_count > min_sup_reads) {
            tr_kept.push_back(tr);
        }
    }

    std::unordered_map<std::string, int> unique_tr_count;
    for (const auto & [read_name, entries] : read_dict) {
        if (entries[0].tr_cov > 0.9) {
            unique_tr_count[entries[0].tr]++;
        }
    }

    for (const auto & [read_name, entries] : read_dict) {
        // collect the entries from the read_dict that have enough coverage to have been retained
        std::vector<ReadDictEntry> tmp;
        for (const auto & entry : entries) {
            if (std::count(tr_kept.begin(), tr_kept.end(), entry.tr) > 0) {
                tmp.push_back(entry);
            }
        }

        ReadDictEntry hit;
        if (tmp.size() > 0) {
            hit = tmp[0];
        } else { // if there are none left, just update 2 and skip it
            count_stat[2] += 1;
            continue;
        }

		// shoud we implement safety checking on this split?
        std::pair<std::string, std::string> bc_umi_split = 
			parseDelim(parseDelim(read_name, '#').first, '_');

        std::string bc = bc_umi_split.first;
        std::string umi = bc_umi_split.second;

        if (kwargs.count("bc_file")) {
            bc = bc_dict[bc];
        }

        if (tmp.size() == 1 && tmp[0].quality > 0) {
            // and then add the umi
            bc_tr_count_dict[bc][hit.tr].push_back(umi);
            count_stat[3] += 1;
        } else if (tmp.size() > 1 &&
            tmp[0].AS_tag == tmp[1].AS_tag &&
            tmp[0].length == tmp[1].length) {
            if (hit.AS_tag > 0.8) {
                bc_tr_count_dict[bc][hit.tr].push_back(umi);
                count_stat[3] += 1;
            } else {
				bc_tr_badcov_count_dict[bc][hit.tr].push_back(umi);
                count_stat[4] += 1;
            }
        } else if (hit.tr_cov < min_tr_coverage ||
            hit.length < min_read_coverage) {
            bc_tr_badcov_count_dict[bc][hit.tr].push_back(umi);
			count_stat[1] += 1;
        } else {
            bc_tr_count_dict[bc][hit.tr].push_back(umi);
            count_stat[3] += 1;
        }
    }
    
	Rcpp::Rcout 
		<< "\tMapping counts:\n"
		<< "\t\tUnmapped: " << count_stat[0]
		<< "\n\t\tNot enough coverage: " << count_stat[1]
		<< "\n\t\tNo good match: " << count_stat[2]
		<< "\n\t\tCounted reads: " << count_stat[3]
		<< "\n\t\tAmbiguous reads: " << count_stat[4]
		<< "\n\t\tNot in annotation: " << count_stat[5] << "\n";
    return RealignedBamData {bc_tr_count_dict, bc_tr_badcov_count_dict, tr_kept};
}

/*
    a debugging function,
    logs the contents of realignedBamData to make sure eveything is being populated as intended
*/
void
log_realigned(RealignedBamData realignedBamData)
{
    std::cout << "started log_realigned\n";

    std::cout << "bc_tr_badcov_count_dict (size " << realignedBamData.bc_tr_badcov_count_dict.size() << "):\n";
    for (const auto & [key, val] : realignedBamData.bc_tr_badcov_count_dict) {
        std::cout << "\t" << key << ": (size " << val.size() << ") [";
        for (const auto & v : val) {
            std::cout << "(" << v.first << "," << v.second.size() << ") ";
        }
        std::cout << "]\n";
    }
    std::cout << "\nbc_tr_count_dict (size " << realignedBamData.bc_tr_count_dict.size() << "):\n";
    for (const auto & [key, val] : realignedBamData.bc_tr_count_dict) {
        std::cout << "\t" << key << ": (size " << val.size() << ") [";
        for (const auto & v : val) {
            std::cout << "(" << v.first << "," << v.second.size() << ") ";
        }
        std::cout << "]\n";
    }
    std::cout << "\ntr_kept (size " << realignedBamData.tr_kept.size() << "):\n";
    for (const auto & k : realignedBamData.tr_kept) {
        std::cout << "\t" << k << "\n";
    }
}


std::unordered_map<std::string, std::string>
make_bc_dict(std::string bc_anno)
{
    /*  
        opens a bc file, parses it into a dictionary format
    */

    // open the file
    std::ifstream
    bc_file (bc_anno);

    // make a dictionary to output to
    std::unordered_map<std::string, std::string>
    bc_dict;

    // now parse the file line by line
    std::string line;
    bool header_line = true;
    while (std::getline(bc_file, line)) {
        // skip the first line
        if (header_line) {
            header_line = false;
            continue;
        }

        std::stringstream
        linestream (line);

        // break each line into words
        std::vector<std::string> words;
        std::string word;
        while (std::getline(linestream, word, ' ')) {
            words.push_back(word);
        }

        bc_dict[words[1]] = words[0];
    }

    return bc_dict;
}