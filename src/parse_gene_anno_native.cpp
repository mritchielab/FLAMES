#include "parse_gene_anno_native.h"

using namespace Rcpp;

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

    // for every record in the GFF3 file
    GFFRecord rec = parser.nextRecord();
    int recnum = -5;
    while (!parser.empty()) {
        std::cout << "\tparsing a record\n";
        for (const auto & [key, val] : rec.attributes) {
            std::cout << "\t\t" << key << ":" << val << "\n";
        }
        // std::cout << "\tparsing record " << recnum << "\n";
        recnum++;

        if (rec.type == "gene") {
            // add this records data to chr_to_gene
            gff_data.chr_to_gene[rec.seqid].push_back(rec.attributes["gene_id"]);
        } else if (rec.type == "transcript") {
            if (rec.attributes["gene_id"].length() == 0) {
                // std::cout << "\t\t\trec has no gene_id\n";
                warning("Transcript did not have 'gene_id' attribute: %s", parser.formatGFFRecordAttributes(rec));
                continue;
            }

            std::string gene_id = rec.attributes["gene_id"];
            // insert the gene_id into the chr_to_gene[rec.seqid] map.
            // the chosen DS ensures there are no duplicates.
            if (std::count(gff_data.chr_to_gene[rec.seqid].begin(), gff_data.chr_to_gene[rec.seqid].end(), gene_id) == 0) {
                // std::cout << "\t\t\tappending " << gene_id << " to chr_to_gene[" << rec.seqid << "]\n";
                gff_data.chr_to_gene[rec.seqid].push_back(gene_id);
            }
            // std::cout << "\t\t\tappending " << rec.attributes["transcript_id"] << " to gene_to_transcript[" << gene_id << "]\n";
            gff_data.gene_to_transcript[gene_id].push_back(rec.attributes["transcript_id"]);
            // std::cout << "\t\t\tappending {" << rec.start-1 << "," << rec.end << "} to transcript_dict[" << rec.attributes["transcript_id"] << "]\n";
            gff_data.transcript_dict[rec.attributes["transcript_id"]] = {rec.seqid, rec.start - 1, rec.end, rec.strand[0], gene_id};

        } else if (rec.type == "exon") {
            // std::cout << "\t\trec type was exon\n";
            if (rec.attributes["gene_id"].length() == 0) {
                // std::cout << "\t\t\trec has no gene_id\n";
                warning("Exon did not have 'gene_id' attribute: %s", parser.formatGFFRecordAttributes(rec));
                continue;
            }

            if (std::count(gff_data.chr_to_gene[rec.seqid].begin(), gff_data.chr_to_gene[rec.seqid].end(), rec.attributes["gene_id"]) == 0) {
                // std::cout << "\t\t\tappending " << rec.attributes["gene_id"] << " to chr_to_gene[" << rec.seqid << "]\n";
                gff_data.chr_to_gene[rec.seqid].push_back(rec.attributes["gene_id"]);
            }
            
            // if "transcript_id" in rec.attributes keys
            if (rec.attributes["transcript_id"].length() != 0) {
                // std::cout << "\t\t\ttranscript_id in rec.attributes\n";

                // if rec.attributes["transcript_id"] not in transcript_dict
                if (gff_data.transcript_dict[rec.attributes["transcript_id"]].chr.length() == 0) {
                    // std::cout << "\t\t\tappending " << rec.attributes["transcript_id"] << " to gene_to_transcript[" << rec.attributes["gene_id"] << "]\n";
                    gff_data.gene_to_transcript[rec.attributes["gene_id"]].push_back(rec.attributes["transcript_id"]);
                    // std::cout << "\t\t\tappending (" << -1 << "," << 1 << ") to transcript_dict[" << rec.attributes["gene_id"] << "]\n";
                    gff_data.transcript_dict[rec.attributes["transcript_id"]] = {rec.seqid, -1, 1, rec.strand[0], rec.attributes["gene_id"]};
                }

                StartEndPair sep {rec.start - 1, rec.end};
                // std::cout << "\t\t\tappending (" << sep.start << "," << sep.end << ") to transcript_to_exon[" << rec.attributes["transcript_id"] << "]\n";
                gff_data.transcript_to_exon[rec.attributes["transcript_id"]].push_back(sep);
            } else {
                warning("Exon did not have 'transcript_id' attribute: %s", parser.formatGFFRecordAttributes(rec));
            }
        }

        // after processing is done, get the next record
        rec = parser.nextRecord();
    }

    std::ofstream
    outfile ("transcript_dict.txt");
    outfile << "\ntranscript_dict before removing dupes:\n";
    for (const auto & [tr, pos] : gff_data.transcript_dict) {
        outfile << "\t" << tr << ":\n\t\t(" << pos.start << "," << pos.end << ")\n";
    }
    // Remove duplicates from the transcript_to_exon map
    gff_data.remove_transcript_duplicates(true);
    outfile << "\ntranscript_dict after removing dupes:\n";
    for (const auto & [tr, pos] : gff_data.transcript_dict) {
        outfile << "\t" << tr << ":\n\t\t(" << pos.start << "," << pos.end << ")\n";
    }
    parser.close();
    
    std::cout << "finished parse_gtf_tree\n";
    return gff_data;
}


GFFData
parse_gff_tree(std::string filename)
{
    std::cout << "started parse_gff_tree on " << filename << "\n";
    // create an object to store chr_to_gene, transcript_dict, 
    // gene_to_transcript, transcript_to_exon 
    GFFData gff_data;

    std::string annotation_source = guess_annotation_source(filename.c_str());
    std::cout << annotation_source << "\n";

    if (annotation_source == "Ensembl") {
        // create the GFF3 parser
        ParseGFF3 parser (filename.c_str());
        // for every record in the GFF3 file.
        std::cout << "parsing\n";
        GFFRecord rec = parser.nextRecord(true);
        while (!parser.empty()) {
            // std::cout << "\tparsing a record\n";
            for (const auto & [key, val] : rec.attributes) {
                // std::cout << "\t\t" << key << ":" << val << "\n";
            }

            if (rec.attributes.count("gene_id") > 0) {
                gff_data.chr_to_gene[rec.seqid].push_back(rec.attributes["gene_id"]);
                // std::cout << "adding " << rec.attributes["gene_id"] << " to chr_to_gene\n";
            }
            
            if (rec.type == "transcript") {
                std::string parent_att = rec.attributes["Parent"];
                int find_colon = parent_att.find(":");

                std::string split0 = parent_att.substr(0, find_colon);
                std::string gene_id = parent_att.substr(find_colon + 1, parent_att.length());

                gff_data.gene_to_transcript[gene_id].push_back(rec.attributes["transcript_id"]);
                Pos pos = {rec.seqid, rec.start-1, rec.end, rec.strand[0], gene_id};
                // std::cout << "adding " << pos.start <<","<< pos.end <<" to transcript_dict\n";
                gff_data.transcript_dict[rec.attributes["transcript_id"]] = pos;
            } else if (rec.type == "exon") {
                // std::cout << "\t\texon\n";
                std::string parent_att = rec.attributes["Parent"];
                int find_colon = parent_att.find(":");

                std::string split0 = parent_att.substr(0, find_colon);
                std::string gene_id = parent_att.substr(find_colon + 1, parent_att.length());

                if (split0 != "transcript") {
                    warning("Format Error.");
                }

                StartEndPair sep {rec.start-1, rec.end};
                gff_data.transcript_to_exon[gene_id].push_back(sep);
                // std::cout << "Ensembl added " << sep.start << "," << sep.end << " to " << gene_id << "\n";
            }

            rec = parser.nextRecord(true);
        }

        parser.close();
    } else if (annotation_source == "GENCODE") {
        // create the GFF3 parser
        ParseGFF3 parser (filename.c_str());
        // for every record in the GFF3 file.
        GFFRecord rec = parser.nextRecord();

        while (!parser.empty()) {
            // std::cout << "rec parent is " << rec.attributes["Parent"] << "\n"; 
            
            if (rec.type == "gene") {
                gff_data.chr_to_gene[rec.seqid].push_back(rec.attributes["gene_id"]); // add the gene_id to the chr_to_gene[seqid] map
            } else if (rec.type == "transcript") {
                std::string gene_id = rec.attributes["gene_id"];

                gff_data.gene_to_transcript[rec.seqid].push_back(gene_id);
                // std::cout << "adding " << rec.start-1 <<","<< rec.end <<" to transcript_dict\n";
                gff_data.transcript_dict[rec.attributes["transcript_id"]] = {rec.seqid, rec.start-1, rec.end, rec.strand[0], gene_id};
            } else if (rec.type == "exon") {
                StartEndPair sep {rec.start-1, rec.end};
                gff_data.transcript_to_exon[rec.attributes["gene_id"]].push_back(sep);
                // std::cout << "GENCODE added " << sep.start << "," << sep.end << " to " << rec.attributes["gene_id"] << "\n";
            }

            rec = parser.nextRecord();
        }

        std::cout << "made it to the end of the parser\n";
        parser.close();
    }

    //gff_data.remove_transcript_duplicates(true);

    // std::cout << "finished parse_gff_tree\n";
    // std::cout << "transcript_dict.size():" << gff_data.transcript_dict.size() << "\n";
    return gff_data;
}

std::unordered_map<std::string, std::vector<StartEndPair>>
GFFData::remove_transcript_duplicates(bool update_transcript_dict)
{
    std::cout << "started remove_transcript_duplicates\n";
    std::unordered_map<std::string, std::vector<StartEndPair>>
    transcript_to_exon_new;

    // Remove duplicates from the transcript_to_exon maps
    for (auto tr = this->transcript_to_exon.begin(); tr != this->transcript_to_exon.end(); ++tr) {
        // it->first is std::string key
        // it->second is std::vector<StartEndPair> list of pairs
        std::sort(tr->second.begin(), tr->second.end(), StartEndPairCompare);
        // if the forward list has at least two elements and start of first element equals start of second
        if (tr->second.size() >= 2 &&
            (tr->second.begin()->start) == ((tr->second.begin())+1)->start) {
            
            std::vector<StartEndPair>
            pair_vec = {};
            pair_vec.push_back(tr->second.front());
            
            // iterate over every exon in the map, removing duplicates
            // for (const auto & pair : pair_list) {
            //     if (pair.start != (int)(++(pair_list.end()))[0] and pair.end != (int)(++(pair_list.begin()))[0]) {
            //         transcript_to_exon_new[tr->first].insert(pair);
            //     }
            // }

            for (auto ex = tr->second.begin() + 1; ex != tr->second.end(); ++ex) {
                if (ex->start != (pair_vec.end())->start && ex->end != (pair_vec.begin())->start) {
                    if (transcript_to_exon_new.count(tr->first) == 0) {
                        transcript_to_exon_new[tr->first] = {};
                    }
                    transcript_to_exon_new[tr->first].push_back({ex->start, ex->end});
                }
            }
            transcript_to_exon_new[tr->first] = pair_vec;

            if (update_transcript_dict) {
                // update transcript_dict[this transcript] with a new Pos object of
                // the correct start and end positions of this exon
                this->transcript_dict[tr->first] = {
                    this->transcript_dict[tr->first].chr,
                    pair_vec.begin()->start,
                    pair_vec.end()->end,
                    this->transcript_dict[tr->first].strand,
                    this->transcript_dict[tr->first].parent_id
                };
            }
        } else {
            // If there aren't enough elements to remove duplicates (more than 2),
            // Just simply add them to the new list and update transcript_dict

            std::vector<StartEndPair>
            pair_vec = {};

            for (auto ex = tr->second.begin(); ex != tr->second.end(); ++ex) {
                pair_vec.push_back({ex->start, ex->end});
            }

            transcript_to_exon_new[tr->first] = pair_vec;

            if (update_transcript_dict) {
                this->transcript_dict[tr->first] = Pos {
                    this->transcript_dict[tr->first].chr,
                    this->transcript_to_exon[tr->first].begin()->start,
                    (this->transcript_to_exon[tr->first].end()-1)->end,
                    this->transcript_dict[tr->first].strand,
                    this->transcript_dict[tr->first].parent_id
                };
            }
        }
    }

    std::cout << "finished remove_transcript_duplicates\n";
    return transcript_to_exon_new;
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

void
log_gff_data
(
    GFFData gff_data, 
    std::string filename
)
{
    /*
        records all the details of gff_data in a file
    */

    std::ofstream
    file (filename);

    file << "chr_to_gene: (size " << gff_data.chr_to_gene.size() << ")\n";
    for (const auto & [chr, gene] : gff_data.chr_to_gene) {
        file << "\t" << chr << ": (size " << gene.size() << ") [";
        for (const auto & g : gene) {
            file << g << ", ";
        }
        file << "]\n";
    }

    file << "\ntranscript_dict: (size " << gff_data.transcript_dict.size() << ")\n";
    for (const auto & [tr, pos] : gff_data.transcript_dict) {
        file << "\t" << tr << ":(" 
            << pos.chr << "," 
            << pos.start << "," 
            << pos.end << "," 
            << pos.parent_id << "," 
            << pos.strand << ")\n";
    }

    file << "\ngene_to_transcript: (size " << gff_data.gene_to_transcript.size() << ")\n";
    for (const auto & [gene, transcript] : gff_data.gene_to_transcript) {
        file << "\t" << gene << ": (size " << transcript.size() << ") [";
        for (const auto & tr : transcript) {
            file << tr << ", ";
        }
        file << "]\n";
    }

    file << "\ntranscript_to_exon: (size " << gff_data.transcript_to_exon.size() << ")\n";
    for (const auto & [transcript, exon] : gff_data.transcript_to_exon) {
        file << "\t" << transcript << ": (size " << exon.size() << ") [";
        for (const auto & ex : exon) {
            file << "(" << ex.start << "," << ex.end << "), ";
        }
        file << "]\n";
    }
}