#include "gff3_to_fa.h"

/*  read data from fa_file,
    check for duplicates,
    log the newly updated content to fa_out_f
*/
void
get_transcript_seq
(
    std::string fa_file,
    std::string fa_out_f,
    std::unordered_map<std::string, std::vector<std::string>>   * chr_to_gene,
    std::unordered_map<std::string, Pos>                        * transcript_dict,
    std::unordered_map<std::string, std::vector<std::string>>   * gene_to_transcript,
    std::unordered_map<std::string, std::vector<StartEndPair>>  * transcript_to_exon,
    ReferenceDict * ref_dict
)
{
    std::cout << "started get_transcript_seq\n";

    std::unordered_map<std::vector<StartEndPair>, std::string>
    global_isoform_dict;
    std::unordered_map<std::string, std::string>
    global_seq_dict = {};
    std::unordered_map<std::string, std::string>
    fa_dict = {};

    // open the output file
    std::ofstream 
    fa_out(fa_out_f);
    
    // load in data from the FASTA input
    std::unordered_map<std::string, std::string>
    raw_dict = get_fa_simple(fa_file);

    // then look through all the data we just loaded in
    std::cout << "iterating raw_dict, size " << raw_dict.size() << "\n";
    for (const auto & [chr, seq] : raw_dict) {
        std::cout << "current iter:" << chr <<","<<seq<<"\n";
        // first, check that the chr is in chr_to_gene
        if ((*chr_to_gene).find(chr) == (*chr_to_gene).end()) {
            std::cout << "skipping\n";
            continue;
        }

        for (const auto & gene : (*chr_to_gene)[chr]) {
            for (const auto & transcript : (*gene_to_transcript)[gene]) {
                // make a list of every StartEndPair in this transcript
                std::vector<StartEndPair>
                iso_list = {};
                for (const auto & exon : (*transcript_to_exon)[transcript]) {
                    iso_list.push_back(exon);
                }

                // check that this exact StartEndPair list isn't already in global_isoform_dict
                if (global_isoform_dict.find(iso_list) != global_isoform_dict.end()) {
                    std::cout << "Duplicate transcript annotation: " << global_isoform_dict[iso_list] << ", " << transcript << "\n";
                } else {
                    global_isoform_dict[iso_list] = transcript;

                    // now, build a string that is the full sequence of the given transcript
                    std::string
                    transcript_seq = "";
                    for (const auto & exon : (*transcript_to_exon)[transcript]) {
                        transcript_seq.append(seq.substr(exon.start, exon.end));
                    }

                    // check if we need to reverse and swap the strand
                    if ((*transcript_dict)[transcript].strand != '+') {
                        transcript_seq = r_c(&transcript_seq);
                    }

                    fa_dict[transcript] = transcript_seq;
                    if (global_seq_dict.find(transcript_seq) != global_seq_dict.end()) {
                        std::cout << "Duplicate transcript sequence: " << global_seq_dict[transcript_seq] << ", " << transcript << "\n";
                    } else {
                        global_seq_dict[transcript_seq] = transcript;
                    }
                }
            }

            if (ref_dict != nullptr) {
                // if the reference dictionary contains the chr
                if (ref_dict->chr_to_gene.find(chr) != ref_dict->chr_to_gene.end()) {
                    for (const auto & transcript : (*ref_dict).gene_to_transcript[gene]) {
                        // first check that the transcript isn't in transcript_to_exon
                        if ((*transcript_to_exon).find(transcript) != (*transcript_to_exon).end()) {
                            continue;
                        }

                        std::vector<StartEndPair>
                        iso_list = {};

                        for (const auto & exon : ref_dict->transcript_to_exon[transcript]) {
                            iso_list.push_back(exon);
                        }

                        // check to see if the transcript is in global_isoform_dict
                        if (global_isoform_dict.find(iso_list) != global_isoform_dict.end()) {
                            if (ref_dict->transcript_to_exon.find(global_isoform_dict[iso_list]) != ref_dict->transcript_to_exon.end()) {
                                std::cout << "Transcript with the same coordination: " << global_isoform_dict[iso_list] << ", " << transcript << "\n";
                                global_seq_dict[fa_dict[global_isoform_dict[iso_list]]] = transcript;
                            }
                        } else {
                            global_isoform_dict[iso_list] = transcript;
                            
                            std::string
                            transcript_seq;
                            for (const auto & exon : ref_dict->transcript_to_exon[transcript]) {
                                transcript_seq.append(seq.substr(exon.start, exon.end));
                            }

                            // reverse and switch strands if we need to
                            if (ref_dict->transcript_dict[transcript].strand != '+') {
                                transcript_seq = r_c(&transcript_seq);
                            }
                            
                            // log a duplicate if we need to
                            if (global_seq_dict.find(transcript_seq) != global_seq_dict.end()) {
                                std::cout << "Duplicate transcript sequence: " << global_seq_dict[transcript_seq] << ", " << transcript << "\n";
                            }
                            
                            global_seq_dict[transcript_seq] = transcript;
                        }
                    }
                }
            }
        }
    }

    for (const auto & [transcript_seq, transcript] : global_seq_dict) {
        write_fa(&fa_out, transcript, transcript_seq);
        std::cout << "written line to fa_out\n";
    }
    fa_out.close();
}

/*  take a FASTA file, 
    parse it and return a map of chr header names to full sequence strings
*/
std::unordered_map<std::string, std::string>
get_fa_simple
(
    std::string filename
)
{
    /* a quick lambda function */
    auto first_space = [] (std::string line) {
        for (int i = 0; i < line.length(); ++i) {
            if (line[i] == ' ') {
                return i;
            }
        }
        return (int)line.length();
    };
    std::unordered_map<std::string, std::string>
    output = {};

    std::string
    chr = "";
    std::string
    full_seq = "";

    std::ifstream
    fa_in(filename);

    std::string line;
    while (std::getline(fa_in, line)) {

        if (line[0] == '>') { // a new sequence has started
            // so push the old one
            if (chr != "") {
                output[chr] = full_seq;
            }
            
            // and start a new sequence
            chr = line.substr(1,first_space(line) - 1);
            full_seq = "";
        } else { // just a regular line
            full_seq.append(line);
        }
    }
    
    if (output.count(chr) == 0) {
        output[chr] = full_seq;
    }

    fa_in.close();
    return output;
}

/*  takes an opened fa_out file, and writes a new entry to it character by character,
    wrapping after each wrap_len characters
*/
void
write_fa(std::ofstream * fa_out, std::string na, std::string seq, int wrap_len)
{
    (*fa_out) << ">" << na << "\n";

    for (int i = 0; i < seq.length(); ++i) {
        // break the line if we need to
        if ((i > 0) && (i % wrap_len == 0)) {
            (*fa_out) << "\n";
        }

        // log the next character
        (*fa_out) << seq[i];
    }
    (*fa_out) << "\n";
}

/*  takes a full sequence,
    reverses it and swaps all of the characters using CP
*/
std::string
r_c(std::string * seq)
{
    std::unordered_map<char, char>
    CP = {
        {'A', 'T'},
        {'T', 'A'},
        {'C', 'G'},
        {'G', 'C'},
        {'N', 'N'},
        {'a', 't'},
        {'t', 'a'},
        {'c', 'g'},
        {'g', 'c'}
    };

    std::string
    new_seq = "";

    for (int i = (int)(*seq).length() - 1; i >= 0; --i) {
        new_seq.push_back(CP[(*seq)[i]]);
    }
    return new_seq;
}