#include "get_transcript_seq.h"

#include <string>
#include <sstream>
#include <unordered_map>
#include <vector>
#include <fstream>
#include <Rcpp.h>
#include "zlib.h"

#include "../classes/GFFData.h"


/*  take a FASTA file, 
    parse it and return a map of chr header names to full sequence strings
*/
std::unordered_map<std::string, std::string>
get_fa(const std::string &filename)
{
    /* a quick lambda function to find the position of the first space in a line */
    auto first_space = [] (const std::string &line) {
        for (int i = 0; i < (int)line.length(); ++i) {
            if (line[i] == ' ') {
                return i;
            }
        }
        return (int)line.length();
    };

    std::unordered_map<std::string, std::string> output;

    std::string chr = "";
    std::string full_seq = "";

    const int MAXLEN = 65535;
    char buf[MAXLEN];
    gzFile fa_in = gzopen(filename.c_str(), "r");

    std::string line;
    while (gzgets(fa_in, buf, MAXLEN)) {
        line = std::string(buf);

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
    // push the last one
    if (chr != "") {
        output[chr] = full_seq;
    }

    gzclose(fa_in);
    return output;
}

std::string
reverse_complement(const std::string &seq) 
{
    const std::unordered_map<char, char> CP {
        {'A', 'T'}, {'T', 'A'}, {'C', 'G'}, {'G', 'C'}, {'N', 'N'},
        {'a', 't'}, {'t', 'a'}, {'c', 'g'}, {'g', 'c'}
    };
    std::stringstream ss;
    for (int i = (int)seq.size(); i >= 0; i--) {
        ss << CP.at(seq[i]);
    }
    return ss.str();
}

void
write_fa(const std::string &fa_out_file, const std::string &na, const std::string &seq, int wrap_len = 50)
{
    // open the output file
    std::ofstream fa_out(fa_out_file);

    fa_out << ">" << na << "\n";
    int i = 0;
    while (i < seq.size()) {
        if (i + wrap_len > seq.size()) {
            fa_out << seq.substr(i, seq.size() - i) << "\n";
            break;
        }

        fa_out << seq.substr(i, wrap_len) << "\n";
        i += wrap_len;
    }

    fa_out.close();
}


void
get_transcript_seq(
    const std::string &fa_file,
    const std::string &fa_out_file,
    const GFFData &isoform_annotation,
    const GFFData &ref_annotation)
{
    std::unordered_map<std::vector<int>, std::string>   global_isoform_dict;
    std::unordered_map<std::string, std::string>        global_seq_dict;
    std::unordered_map<std::string, std::string>        fa_dict;
    
    // load in data from the FASTA input
    std::unordered_map<std::string, std::string> raw_dict = get_fa(fa_file);
    
    // then look through all the data we just loaded in
    for (const auto & [chr, seq] : raw_dict) {
        // first, check that the chr is in chr_to_gene
        if (isoform_annotation.chr_to_gene.find(chr) == isoform_annotation.chr_to_gene.end()) {
            continue;
        }

        for (const auto & gene : isoform_annotation.chr_to_gene.at(chr)) {
            for (const auto & transcript : isoform_annotation.gene_to_transcript.at(gene)) {
                // make a list of every StartEndPair in this transcript
                std::vector<int> iso_list;
                for (const auto & exon : isoform_annotation.transcript_to_exon.at(transcript)) {
                    iso_list.push_back(exon.start);
                    iso_list.push_back(exon.end);
                }

                // check that this exact StartEndPair list isn't already in global_isoform_dict
                if (global_isoform_dict.find(iso_list) != global_isoform_dict.end()) {
                    Rcpp::Rcout << "Duplicate transcript annotation: " << global_isoform_dict[iso_list] << ", " << transcript << "\n";
                } else {
                    global_isoform_dict[iso_list] = transcript;

                    // now, build a string that is the full sequence of the given transcript
                    std::string transcript_seq = "";
                    for (const auto & exon : isoform_annotation.transcript_to_exon.at(transcript)) {
                        transcript_seq.append(seq.substr(exon.start, exon.end - exon.start));
                    }

                    // check if we need to reverse and swap the strand
                    if (isoform_annotation.transcript_dict.at(transcript).strand != '+') {
                        transcript_seq = reverse_complement(transcript_seq);
                    }

                    fa_dict[transcript] = transcript_seq;
                    if (global_seq_dict.find(transcript_seq) != global_seq_dict.end()) {
                        std::cout << "Duplicate transcript sequence: " << global_seq_dict[transcript_seq] << ", " << transcript << "\n";
                    } else {
                        global_seq_dict[transcript_seq] = transcript;
                    }
                }
            }
            if (!ref_annotation.is_empty()) {
                // if the reference dictionary contains the chr
                if (ref_annotation.chr_to_gene.count(chr)) {
                    for (const auto & transcript : ref_annotation.gene_to_transcript.at(gene)) {
                        // first check that the transcript isn't in transcript_to_exon
                        if (isoform_annotation.transcript_to_exon.count(transcript)) {
                            continue;
                        }

                        std::vector<int> iso_list;

                        for (const auto & exon : ref_annotation.transcript_to_exon.at(transcript)) {
                            iso_list.push_back(exon.start);
                            iso_list.push_back(exon.end);
                        }
    
                        // check to see if the transcript is in global_isoform_dict
                        if (global_isoform_dict.find(iso_list) != global_isoform_dict.end()) {
                            if (ref_annotation.transcript_to_exon.count(global_isoform_dict[iso_list])) {
                                Rcpp::Rcout << "Transcript with the same coordination: " << global_isoform_dict[iso_list] << ", " << transcript << "\n";
                                global_seq_dict[fa_dict[global_isoform_dict[iso_list]]] = transcript;
                            }
                        } else {
                            global_isoform_dict[iso_list] = transcript;
                            
                            std::string transcript_seq;
                            for (const auto & exon : ref_annotation.transcript_to_exon.at(transcript)) {
                                transcript_seq.append(seq.substr(exon.start, exon.end - exon.start));
                            }
                            
                            // reverse and switch strands if we need to
                            if (ref_annotation.transcript_dict.at(transcript).strand != '+') {
                                transcript_seq = reverse_complement(transcript_seq);
                            }
                            
                            // log a duplicate if we need to
                            if (global_seq_dict.count(transcript_seq) ) {
                                Rcpp::Rcout << "Duplicate transcript sequence: " << global_seq_dict[transcript_seq] << ", " << transcript << "\n";
                            }
                            
                            global_seq_dict[transcript_seq] = transcript;
                        }
                    }
                }
            }
        }
    }

    // for (const auto & [transcript_seq, transcript] : global_seq_dict) {
    //     write_fa(fa_out_file, transcript, transcript_seq);
    // }
}
