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
    std::cout << "started file_to_map\n";
    std::unordered_map<std::string, int>
    map;

    std::ifstream
    file(filename);

    std::string line;
    int lines = 0;
    while (std::getline(file, line)) {
        std::cout << "line:" <<line<<"\n";
        lines++;
        std::stringstream
        linestream (line);

        // break each line into words
        std::vector<std::string> words;
        std::string word;
        while (std::getline(linestream, word, '\t')) {
            words.push_back(word);
        }
        std::cout << "found " << words.size() << " words\n";
        if (words.size() >= 2) {
            map[words[0]] = atoi(words[1].c_str());
            std::cout << "\t(added to map)\n";
        }
    }

    std::cout << "file was " << lines << " lines long\n";

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
    bam_index_t *bam_index = bam_index_load(bam_in.c_str());
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

int
query_len(std::string cigar_string, bool hard_clipping)
{
    /*
        https://stackoverflow.com/questions/39710796/infer-the-length-of-a-sequence-using-the-cigar
        Given a CIGAR string, return the number of bases consumed from the
        query sequence.
        CIGAR is a sequence of the form <operations><operator> such that operations is an integer giving 
        the number of times the operator is used
        M = match
        I = Insertion
        S = Soft clipping
        = = sequence match
        X = sequence mismatch
    */

   // set up which operations should consume a read
   std::vector<char> read_consuming_ops;
   if (!hard_clipping) {
       read_consuming_ops = {'M', 'I', 'S', '=', 'X'};
   } else {
       read_consuming_ops = {'M', 'I', 'S', 'H', '=', 'X'};
   }

   int result = 0;
   int this_len = 0;
   for (const auto & c : cigar_string) {
       if (isdigit(c)) {
           // if it's an int, update the len
           this_len = int(c);
       } else {
           // check if it's in the read_consuming_ops
           if (std::count(read_consuming_ops.begin(), read_consuming_ops.end(), c) > 0) {
               // if so we have to add it
               result += this_len;
           }
       }
   }

   return result;
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
    std::cout << "started parse_realigned_bam\n";

    // we need to read in the fa_idx_f file line by line, adding each one to the dict
    std::unordered_map<std::string, int>
    fa_idx = file_to_map(fa_idx_f);
    std::cout << "fa_idx is " << fa_idx.size() << " long\n"; 

    std::unordered_map<std::string, std::unordered_map<std::string, std::vector<std::string>>>
    bc_tr_count_dict = {};

    std::unordered_map<std::string, std::unordered_map<std::string, std::vector<std::string>>>
    bc_tr_badcov_count_dict = {};

    std::unordered_map<std::string, std::vector<float>>
    tr_cov_dict = {};

    std::unordered_map<std::string, std::vector<ReadDictEntry>>
    read_dict;
    
    std::unordered_map<std::string, int>
    count_stat = {};
    count_stat["unmapped"] = 0;
    count_stat["not in annotation"] = 0;
    count_stat["no good match"] = 0;
    count_stat["counted_reads"] = 0;
    count_stat["ambiguous_reads"] = 0;
    count_stat["not_enough_coverage"] = 0;


    std::unordered_map<std::string, std::string>
    bc_dict = {};
    if (kwargs.count("bc_file") > 0) {
        bc_dict = make_bc_dict(kwargs["bc_file"]);
    }
    std::cout << "bc_dict is " << bc_dict.size() << " long\n";

    
    // read a bamfile
    bamFile bam = bam_open(bam_in.c_str(), "r"); // bam.h
    bam_index_t *bam_index = bam_index_load(bam_in.c_str());
    bam_header_t *header = bam_header_read(bam); // bam.h

    // produce and populate a records file
    std::vector<BAMRecord>
    records = {};
    // iterate over every entry in the bam
    bam1_t *b = bam_init1();
    while (bam_read1(bam, b) >= 0) {
        if (b->core.tid == -1) {
            continue;
        }
        std::cout << "\t\treading one from bam\n";
        BAMRecord rec = read_record(b, header);
        records.push_back(rec);
    }
    
    bam_close(bam);

    std::cout << "records is " << records.size() << " long\n";
    for (const auto & rec : records) {
        // if it's unmapped, just update the count and continue
        if (rec.flag.read_unmapped) {
            count_stat["unmapped"] += 1;
            continue;
        }

        int
        map_start = rec.reference_start;
        int
        map_end = rec.reference_end;

        std::string
        tr = rec.reference_name;
        
        
        float
        tr_cov = (float)((map_end - map_start)/fa_idx[tr]);
        // if there is no dictionary entry for this tr_cov, create one
        if (tr_cov_dict.count(tr) == 0) {
            tr_cov_dict[tr] = {};
        }
        // then add the new cov
        tr_cov_dict[tr].push_back(tr_cov);

        int inferred_read_length = query_len(rec.cigar_string);
        float tr_length = (float)(rec.query_alignment_length)/inferred_read_length;

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
            read_dict[rec.read_name] = {};
            
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
            count_stat["not in annotation"] += 1;
            std::cout << "\t" << tr << " not in annotation ???\n";
        }
    }

    // build up a new dict of only the tr entries that have sufficient coverage
    // (at least 90% coverage on at least min_sup_reads)
    std::vector<std::string>
    tr_kept;
    for (const auto & [tr, covs] : tr_cov_dict) {
        // count up how many of the coverages are above the acceptable value of 90%
        int sup_read_count = 0;
        for (const auto & cov : covs) {
            if (cov > 0.9) {
                sup_read_count++;
            }
        }
        
        // if the overall threshold is met for this tr, add it to the new dict
        if (sup_read_count > min_sup_reads) {
            tr_kept.push_back(tr);
        }
    }

    std::unordered_map<std::string, int>
    unique_tr_count;
    for (const auto & [read_name, entries] : read_dict) {
        if (entries[0].tr_cov > 0.9) {
            if (unique_tr_count.count(entries[0].tr) == 0) {
                unique_tr_count[entries[0].tr] = 0;
            }
            unique_tr_count[entries[0].tr]++;
        }
    }

    for (const auto & [read_name, entries] : read_dict) {
        // collect the entries from the read_dict that have enough coverage to have been retained
        std::vector<ReadDictEntry>
        tmp;
        for (const auto & entry : entries) {
            if (std::count(tr_kept.begin(), tr_kept.end(), entry.tr) > 0) {
                tmp.push_back(entry);
            }
        }

        ReadDictEntry hit;
        if (tmp.size() > 0) {
            hit = tmp[0];
        } else { // if there are none left, just update "no good match" and skip it
            count_stat["no good match"] += 1;
            continue;
        }

        std::vector<std::string>
        read_name_split;

        std::string current_split = ""; 
        for (const auto & c : read_name) {
            if (c != '#' &&
                c != '_') {
                // just add the character
                current_split += c;
            } else {
                // or add it to the end and reset the current string
                read_name_split.push_back(current_split);
                current_split = "";
            }

            if (c == '#') {
                // if we find one of these, also terminate the loop
                break;
            }
        }

        std::string
        bc = read_name_split[read_name_split.size() - 2];
        std::string
        umi = read_name_split[read_name_split.size() - 1];

        if (kwargs.count("bc_file")) {
            bc = bc_dict[bc];
        }

        if (tmp.size() == 1 &&
            tmp[0].quality > 0) {
            // initialise the dictionary entries if necesary
            if (bc_tr_count_dict.count(bc) == 0) {
                bc_tr_count_dict[bc] = {};
            }
            if (bc_tr_count_dict[bc].count(hit.tr) == 0) {
                bc_tr_count_dict[bc][hit.tr] = {};
            }
            // and then add the umi
            bc_tr_count_dict[bc][hit.tr].push_back(umi);
            count_stat["counted_reads"] += 1;
        } else if (tmp.size() > 1 &&
            tmp[0].AS_tag == tmp[1].AS_tag &&
            tmp[0].tr_cov == tmp[1].tr_cov) {
            if (hit.AS_tag > 0.8) {
                if (bc_tr_count_dict.count(bc) == 0) {
                    bc_tr_count_dict[bc] = {};
                }
                if (bc_tr_count_dict[bc].count(hit.tr) == 0) {
                    bc_tr_count_dict[bc][hit.tr] = {};
                }

                bc_tr_count_dict[bc][hit.tr].push_back(umi);
                count_stat["counted_reads"] += 1;
            } else {
                count_stat["ambiguous_reads"] += 1;
                if (bc_tr_badcov_count_dict.count(bc) == 0) {
                    bc_tr_badcov_count_dict[bc] = {};
                }
                if (bc_tr_badcov_count_dict[bc].count(hit.tr) == 0) {
                    bc_tr_badcov_count_dict[bc][hit.tr] = {};
                }
                bc_tr_badcov_count_dict[bc][hit.tr].push_back(umi);
            }
        } else if (hit.tr_cov < min_tr_coverage ||
            hit.length < min_read_coverage) {
            count_stat["not_enough_coverage"] += 1;
            if (bc_tr_badcov_count_dict.count(bc) == 0) {
                bc_tr_badcov_count_dict[bc] = {};
            }
            if (bc_tr_badcov_count_dict[bc].count(hit.tr) == 0) {
                bc_tr_badcov_count_dict[bc][hit.tr] = {};
            }
            bc_tr_badcov_count_dict[bc][hit.tr].push_back(umi);
        } else {
            if (bc_tr_count_dict.count(bc) == 0) {
                bc_tr_badcov_count_dict[bc] = {};
            }
            if (bc_tr_count_dict[bc].count(hit.tr) == 0) {
                bc_tr_badcov_count_dict[bc][hit.tr] = {};
            }
            bc_tr_count_dict[bc][hit.tr].push_back(umi);
            count_stat["counted_reads"] += 1;
        }
    }
    
    std::cout << "Mapping counts:" << "\n";
    for (const auto & [key, val] : count_stat) {
        std::cout << "\t" << key << ": " << val << "\n";
    }
    return RealignedBamData {bc_tr_count_dict, bc_tr_badcov_count_dict, tr_kept};
}



std::unordered_map<std::string, std::string>
make_bc_dict(std::string bc_anno)
{
    /*  
        opens a bc file, parses it into a dictionary format
    */

    std::cout << "started make_bc_dict\n";

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
        std::cout << "line:" << line << "\n";
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

    std::cout << "finished make_bc_dict\n";
    return bc_dict;
}