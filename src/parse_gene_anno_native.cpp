#include <unordered_map>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <list>
#include <Rcpp.h>

#include "ParseGFF3.hpp"
#include "Parser.h"
#include "Pos.h"
#include "StartEndPair.hpp"

#include "parse_gene_anno_native.h"

using namespace Rcpp;

static bool vectorContains(std::string a, std::vector<std::string> v) {
	for (auto it : v) {
		if (a == it) {
			return true;
		}
	}
	return false;
}

// [[Rcpp::export]]
Rcpp::List
parse_gff_or_gtf_R(std::string filename)
{
    return parse_gff_or_gtf(filename).to_R();
}

GFFData
parse_gff_or_gtf(std::string filename)
{
    if (filename.find(".gtf") != std::string::npos) {
        return parse_gtf_tree(filename);
    } else {
        return parse_gff_tree(filename);
    }
}

GFFData
parse_gtf_tree(std::string filename)
{
    std::cout << "started parse_gtf_tree\n";
    // create an object to store chr_to_gene, transcript_dict, 
    // gene_to_transcript, transcript_to_exon 
    GFFData gff_data;
    // create the GFF3 parser
    ParseGFF3 parser (filename.c_str());

    GFFRecord rec = parser.nextRecord();
    while (!parser.empty()) {
		std::string gene_id = rec.attributes["gene_id"];
        if (rec.type == "gene") {
            std::cout << "\t\trec type was gene\n";
            // add this records data to chr_to_gene
            gff_data.chr_to_gene[rec.seqid].push_back(gene_id);
        } else if (rec.type == "transcript") {
            if (gene_id.length() == 0) {
                warning("Transcript did not have 'gene_id' attribute: %s", parser.formatGFFRecordAttributes(rec));
                continue;
            }

			if (!vectorContains(gene_id, gff_data.chr_to_gene[rec.seqid])) {
            	gff_data.chr_to_gene[rec.seqid].push_back(gene_id);
			}

            gff_data.gene_to_transcript[gene_id].push_back(rec.attributes["transcript_id"]);
            std::cout << "\t\t\tappending {" << rec.start-1 << "," << rec.end << "} to transcript_dict[" << rec.attributes["transcript_id"] << "]\n";
            gff_data.transcript_dict[rec.attributes["transcript_id"]] = {rec.seqid, rec.start - 1, rec.end, rec.strand[0], gene_id};

        } else if (rec.type == "exon") {
            if (gene_id.length() == 0) {
                warning("Exon did not have 'gene_id' attribute: %s", parser.formatGFFRecordAttributes(rec));
                continue;
            } else if (!vectorContains(gene_id, gff_data.chr_to_gene[rec.seqid])) {
            	gff_data.chr_to_gene[rec.seqid].push_back(gene_id);
			}
            
            if (rec.attributes.count("transcript_id")) {
                // if rec.attributes["transcript_id"] not in transcript_dict
                if (!gff_data.transcript_dict.count(rec.attributes["transcript_id"])) {
                    gff_data.gene_to_transcript[gene_id].push_back(rec.attributes["transcript_id"]);
                    gff_data.transcript_dict[rec.attributes["transcript_id"]] = Pos {rec.seqid, -1, 1, rec.strand[0], rec.attributes["gene_id"]};
                }

                gff_data.transcript_to_exon[rec.attributes["transcript_id"]].push_back(StartEndPair {rec.start - 1, rec.end});
            } else {
                warning("Exon did not have 'transcript_id' attribute: %s", parser.formatGFFRecordAttributes(rec));
            }
        }

        // after processing is done, get the next record
        rec = parser.nextRecord();
    }

    
    // Remove duplicates from the transcript_to_exon map
    // gff_data.remove_transcript_duplicates(true);

    parser.close();
    
    std::cout << "finished parse_gtf_tree\n";
    return gff_data;
}


// testing 
// static std::string printAttributes(std::unordered_map<std::string, std::string> a) {
// 	std::stringstream ss;
// 	for (auto it : a) {
// 		ss << it.first << "=" << it.second << ";";
// 	}
// 	return ss.str();
// }
// static std::string printRecord(GFFRecord r) {
// 	std::string attr = printAttributes(r.attributes);

// 	std::stringstream ss;
// 	ss << "GFFRecord(seqid=" << r.seqid
// 		<< ", source=" << r.source 
// 		<< ", type=" << r.type
// 		<< ", start=" << r.start
// 		<< ", end=" << r.end 
// 		<< ", score=" << r.score 
// 		<< ", strand=" << r.strand
// 		<< ", phase=" << r.phase
// 		<< ", attributes={" << attr << "}" << "\n";
// 	return ss.str();
// }

GFFData
parse_gff_tree(std::string filename)
{
    Rcpp::Rcout << "started parse_gff_tree on " << filename << "\n";
    // create an object to store chr_to_gene, transcript_dict, 
    // gene_to_transcript, transcript_to_exon 
    GFFData gff_data;

    std::string annotation_source = guess_annotation_source(filename);
    Rcpp::Rcout << annotation_source << "\n";

    if (annotation_source == "Ensembl") {
        // create the GFF3 parser
        ParseGFF3 parser (filename);
        GFFRecord rec = parser.nextRecord();

        while (!parser.empty()) {
            if (rec.attributes.count("gene_id")) {
                gff_data.chr_to_gene[rec.seqid].push_back(rec.attributes["gene_id"]);
            }

			ParseResult pr = parseKeyValue(rec.attributes["Parent"], ':');
			std::string gene_id = pr.second;

            if (rec.attributes.count("Parent") && pr.first == "gene") {
				gff_data.gene_to_transcript[gene_id].push_back(rec.attributes["transcript_id"]);
				gff_data.transcript_dict[rec.attributes["transcript_id"]] = Pos {rec.seqid, rec.start-1, rec.end, rec.strand[0], gene_id};
				
			} else if (rec.type == "exon") {
                if (pr.first != "transcript") {
					// Rcpp::Rcout << "Format error: " << pr.first << " : " << pr.second << "\n";
                    warning("Format Error.");
                }

                gff_data.transcript_to_exon[gene_id].push_back(StartEndPair {rec.start-1, rec.end});
                // Rcpp::Rcout << "added " << sep.start << "," << sep.end << " to " << gene_id << "\n";
            }

            rec = parser.nextRecord();
        }

        parser.close();
    } else if (annotation_source == "GENCODE") {
        // create the GFF3 parser
        ParseGFF3 parser (filename);
        GFFRecord rec = parser.nextRecord();

        while (!parser.empty()) {
            if (rec.type == "gene") {
                gff_data.chr_to_gene[rec.seqid].push_back(rec.attributes["gene_id"]); // add the gene_id to the chr_to_gene[seqid] map
            } else if (rec.type == "transcript") {
                std::string gene_id = rec.attributes["Parent"];

				if (!vectorContains(gene_id, gff_data.chr_to_gene[rec.seqid])) {
					gff_data.chr_to_gene[rec.seqid].push_back(gene_id);
				}

                gff_data.gene_to_transcript[gene_id].push_back(rec.attributes["transcript_id"]);
                gff_data.transcript_dict[rec.attributes["transcript_id"]] = {rec.seqid, rec.start-1, rec.end, rec.strand[0], gene_id};
            } else if (rec.type == "exon") {
                gff_data.transcript_to_exon[rec.attributes["Parent"]].push_back(StartEndPair {rec.start-1, rec.end});
                // Rcpp::Rcout << "added " << sep.start << "," << sep.end << " to " << rec.attributes["gene_id"] << "\n";
            }

            rec = parser.nextRecord();
        }

        std::cout << "made it to the end of the parser\n";
        parser.close();
    }

    gff_data.remove_transcript_duplicates(false);

    Rcpp::Rcout << "finished parse_gff_tree\n";
    return gff_data;
}


static bool StartEndPairCompare(const StartEndPair &a, const StartEndPair &b) {
    // compare a and b, return true if a is 'less than' b
    // in this case, 'less than' is defined if a.start is less than b.start
    return a.start < b.start;
}

void
GFFData::remove_transcript_duplicates(bool update_transcript_dict) {
    // Remove duplicates from the transcript_to_exon maps
    for (auto tr : transcript_to_exon) {
        // it->first is std::string key
        // it->second is std::vector<StartEndPair> list of pairs
        std::sort(tr.second.begin(), tr.second.end(), StartEndPairCompare);
        // if the forward list has at least two elements and start of first element equals start of second
        if (tr.second.size() >= 2 &&
            (tr.second.begin()->start) == (tr.second.begin() + 1)->start) {
            
            std::vector<StartEndPair> new_ex = {*tr.second.begin()};
			// remove duplicates from the old list of exons stored in transcript_to_exon[x]
            for (auto ex : tr.second) {
				StartEndPair last_ex = *(new_ex.end() - 1);
				if (ex.start != last_ex.start && ex.end != last_ex.start) {
					new_ex.push_back(ex);
				}
            }

            transcript_to_exon[tr.first] = new_ex;
        }

		if (update_transcript_dict) {
                // update transcript_dict[this transcript] with a new Pos object of
                // the correct start and end positions of this exon
                transcript_dict[tr.first] = Pos {
                    transcript_dict[tr.first].chr,
                    tr.second.begin()->start,
                    (tr.second.end() - 1)->end,
                    transcript_dict[tr.first].strand,
                    transcript_dict[tr.first].parent_id
                };
            }
    }

	// remove duplicates from gene_to_transcript
	for (auto ge : gene_to_transcript) {
		std::unordered_map<std::string, int> set;
		std::vector<std::string> new_genes;

		for (auto el : ge.second) {
			set[el] = 1;
		}

		for (auto unique : set) {
			new_genes.push_back(unique.first);
		}

		ge.second = new_genes;
	}
}


List
GFFData::to_R()
{
    // /// Create a Rcpp::List from an unordered_map.
    // /// Specificially used for chr_to_gene and gene_to_transcript in order to export each object
    // /// to the R session for later use.
    // /// Maintains order of data
    // auto create_list_map_to_map = [] (std::unordered_map<std::string, std::unordered_map<std::string, bool>> map)
    // {
    //     List result = List::create();

    //     for (auto it = map.begin; it != map.end(); ++it) {
    //         // iterating over every string to unordered_map in map
    //         CharacterVector inner_list = CharacterVector::create();

    //         // iterating over every string in the internal unordered_map
    //         for (auto in = it->second.begin(); in != it->second.end(); ++in) {
    //             inner_list.push_back(in->first);
    //         }

    //         result.push_back(inner_list, it->first);
    //     }

    //     return result;
    // };

    // auto create_list_map_to_pos = [] (std::unordered_map<std::string, Pos> map)
    // {
    //     List result = List::create();

    //     // result should look like: "transcript" = ["pos"=[]]
    //     for (auto it = map.begin(); it != map.end(); ++it) {
    //         Pos pos = it->second;

    //         result.push_back(pos_to_R(&pos), it->first);
    //     }

    //     return result;
    // };

    // auto create_list_map_to_startendpair = [] (std::unordered_map<std::string, std::list<StartEndPair>> map)
    // {
    //     List transcript_to_exon_list = List::create();

    //     for (auto tr = map.begin(); tr !=map.end(); ++tr) {
    //         // it->first is std::string key
    //         // it->second is std::list<StartEndPair> list of pairs
    //         List pair_list = List::create();

    //         for (auto ex = tr->second.begin(); ex != tr->second.end(); ++ex) {
    //             pair_list.push_back(
    //                 IntegerVector::create(
    //                     ex->start,
    //                     ex->end
    //                 )
    //             );
    //         }

    //         transcript_to_exon_list.push_back(pair_list, tr->first);
    //     }

    //     return transcript_to_exon_list;
    // };

    // return List::create(
    //     _["chr_to_gene"] = create_list_map_to_map(this->chr_to_gene),
    //     _["transcript_dict"] = create_list_map_to_pos(this->transcript_dict),
    //     _["gene_to_transcript"] = create_list_map_to_map(this->gene_to_transcript),
    //     _["transcript_to_exon"] = create_list_map_to_startendpair(this->transcript_to_exon)
    // );
    return List::create();
}

void
GFFData::from_R(List list)
{
    /* import an R list */
}