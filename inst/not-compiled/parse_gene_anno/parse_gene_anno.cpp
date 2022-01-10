#include "parse_gene_anno.h"

using namespace Rcpp;


// THESE FUNCTIONS SUFFER FROM A FEW ERRORS OF LOGIC (totally my fault - Oliver) 
// WHICH HAVE BEEN FIXED IN parse_gene_anno_native.cpp
// Difference is this script yields R objects and _native gives C++ objects


// void test_printing(std::unordered_map<String, std::unordered_map<String, bool>> chr_to_gene,
//     std::unordered_map<String, Pos> transcript_dict,
//     std::unordered_map<String, std::unordered_map<String, bool>> gene_to_transcript,
//     std::unordered_map<String, std::list<StartEndPair>> transcript_to_exon) {
    
//     Rcout << "chr_to_gene:\n";
//     for (auto ch = chr_to_gene.begin(); ch != chr_to_gene.end(); ch++) {
//         Rcout << ch->first.get_cstring() << ":\n";
//         for (auto gene = chr_to_gene[ch->first].begin(); gene != chr_to_gene[ch->first].end(); gene++) {
//             Rcout << "\t" << gene->first.get_cstring() <<"\n";
//         }
//     }

//     Rcout << "transcript_dict:\n";
//     for (auto tr = transcript_dict.begin(); tr != transcript_dict.end(); tr++) {
//         Rcout << tr->first.get_cstring() << ":\n";
//         Pos pos = tr->second;
//         Rcout << "\tPos(chr='" << pos.chr << "\', start=" << pos.start << ", end=" << pos.end << ", strand=\'" << pos.strand << "\', parent_id=\'" << pos.parent_id << "\')\n";
//     }

//     Rcout << "gene_to_transcript:\n";
//     for (auto gene : gene_to_transcript) {
//         Rcout << gene.first.get_cstring() << ":\n";
//         for (auto tr : gene.second) {
//             Rcout << "\t" << tr.first.get_cstring() <<"\n";
//         }
//     }

//     Rcout << "transcript_to_exon:\n";
//     for (auto tr : transcript_to_exon) {
//         Rcout << tr.first.get_cstring() << ":\n";
//         for (auto pair : tr.second) {
//             Rcout << pair.start << " " << pair.end << "\n";
//         }
//     }
// }

/// Create a Rcpp::List from an unordered_map.
/// Specificially used for chr_to_gene and gene_to_transcript in order to export each object
/// to the R session for later use.
/// Maintains order of data
List create_list_map_to_map(std::unordered_map<String, std::unordered_map<String, bool>> &map) {
    List result = List::create();

    for (auto it = map.begin(); it != map.end(); ++it) {
        // iterating over every String to unordered_map in map
        CharacterVector inner_list = CharacterVector::create();

        // iterating over every String in the internal unordered_map
        for (auto in = it->second.begin(); in != it->second.end(); ++in) {
            inner_list.push_back(in->first);
        }

        result.push_back(inner_list, it->first);
    }

    return result;
}

/// Create a Rcpp::List from an unordered_map of Pos structs
/// Specifically used for transcript_dict to export to R session
/// Unfortunately, does not give R session a Pos object, instead Pos are rendered as
/// named lists of elements
List create_list_map_to_pos(std::unordered_map<String, Pos> &map) {
    List result = List::create();

    // result should look like: "transcript" = ["pos"=[]]
    for (auto it = map.begin(); it != map.end(); ++it) {
        Pos pos = it->second;

        result.push_back(
            List::create(_["Pos"] = 
                List::create(
                    _["chr"] = pos.chr,
                    _["start"] = pos.start,
                    _["end"] = pos.end,
                    _["strand"] = pos.strand,
                    _["parent_id"] = pos.parent_id)
            ), 
            it->first);
    }

    return result;
}

/// Remove duplicates from the transcript_to_exon map and create a Rcpp::List for export
/// back to the R session.
/// Only removes duplicates if more than 2 elements in the map
List remove_transcript_duplicates_to_list(std::unordered_map<String, std::list<StartEndPair>> &transcript_to_exon, std::unordered_map<String, Pos> &transcript_dict, bool update_transcript_dict) {
    List transcript_to_exon_list = List::create();

    // Remove duplicates from the transcript_to_exon maps and create the Rcpp::List output object
    for (auto tr = transcript_to_exon.begin(); tr != transcript_to_exon.end(); ++tr) {
        // it->first is String key
        // it->second is std::list<StartEndPair> list of pairs
        tr->second.sort(StartEndPairCompare);
        // if the forward_list has at least two elements and start of first element equals start of second
        if (tr->second.size() >= 2 && 
            (tr->second.begin()->start) == (++(tr->second.begin()))->start) {
            
            // pairList holds the new list after removing duplicates
            List pairList = List::create(
                IntegerVector::create(
                    tr->second.begin()->start,
                    tr->second.begin()->end
                )
            );
            
            // iterate over every exon in the map, removing duplicates and adding
            // each pair to the pairList list
            for (auto ex = ++(tr->second.begin()); ex != tr->second.end(); ++ex) {
                if (ex->start != (int)(++(pairList.end()))[0] and ex->end != (int)(++(pairList.begin()))[0]) {
                    pairList.push_back(
                        IntegerVector::create(
                            ex->start,
                            ex->end
                        )
                    );
                }
            }
            transcript_to_exon_list.push_back(pairList, tr->first);

            if (update_transcript_dict) {
                // update transcript_dict[this transcript] with a new Pos object of
                // the correct start and end positions of this exon
                transcript_dict[tr->first] = Pos {
                    transcript_dict[tr->first].chr,
                    pairList.begin()[0],
                    (--(pairList.end()))[1],
                    transcript_dict[tr->first].strand,
                    transcript_dict[tr->first].parent_id};
            }

        } else {
            // If there aren't enough elements to remove duplicates (more than 2),
            // Just simply add them to the new Rcpp::list and update transcript_dict
            List pairList = List::create();

            for (auto ex = tr->second.begin(); ex != tr->second.end(); ++ex) {
                pairList.push_back(
                    IntegerVector::create(
                        ex->start,
                        ex->end
                    )
                );
            }

            transcript_to_exon_list.push_back(pairList, tr->first);

            if (update_transcript_dict) {
                transcript_dict[tr->first] = Pos {
                    transcript_dict[tr->first].chr, 
                    transcript_to_exon[tr->first].begin()->start,
                    (--transcript_to_exon[tr->first].end())->end,
                    transcript_dict[tr->first].strand,
                    transcript_dict[tr->first].parent_id};
            }
        }
    }

    return transcript_to_exon_list;
}

List parse_gff_tree(const char * gff_filename) {
    std::cout << "started parse_gff_tree\n";
        std::unordered_map<String, std::unordered_map<String, bool>>    chr_to_gene; // map of chr to map of genes (second map for efficient searching)
    std::unordered_map<String, Pos>                                 transcript_dict;
    std::unordered_map<String, std::unordered_map<String, bool>>    gene_to_transcript;
    std::unordered_map<String, std::list<StartEndPair>>     transcript_to_exon;

    std::string annotation_source = guess_annotation_source(gff_filename);

    if (annotation_source == "Ensembl") {
        // create the GFF3 parser
        ParseGFF3 parser (gff_filename);
        // for every record in the GFF3 file.

        GFFRecord rec = parser.nextRecord();
        while (!parser.empty()) {
            if (rec.attributes["gene_id"].length() != 0) {
                chr_to_gene[rec.seqid][rec.attributes["gene_id"]] = true;
            }
            if (rec.attributes["Parent"].length() != 0) {
                std::string parent_att = rec.attributes["Parent"];
                int find_colon = parent_att.find(":");

                std::string split0 = parent_att.substr(0, find_colon);
                std::string gene_id = parent_att.substr(find_colon + 1, parent_att.length());

                gene_to_transcript[gene_id][rec.attributes["transcript_id"]] = true;
                Pos pos {rec.seqid, rec.start-1, rec.end, rec.strand[0], gene_id};
                transcript_dict[rec.attributes["transcript_id"]] = pos;
            } else if (rec.type == "exon") {
                std::string parent_att = rec.attributes["Parent"];
                int find_colon = parent_att.find(":");

                std::string split0 = parent_att.substr(0, find_colon);
                std::string gene_id = parent_att.substr(find_colon + 1, parent_att.length());

                if (split0 != "transcript") {
                    warning("Format Error.");
                }
                
                StartEndPair sep {rec.start-1, rec.end};
                transcript_to_exon[gene_id].push_back(sep);
            }

            rec = parser.nextRecord();
        }

        parser.close();

    } else if (annotation_source == "GENCODE") {
        // create the GFF3 parser
        ParseGFF3 parser (gff_filename);
        // for every record in the GFF3 file.
        GFFRecord rec = parser.nextRecord();

        while (!parser.empty()) {
            if (rec.type == "gene") {
                chr_to_gene[rec.seqid][rec.attributes["gene_id"]] = true; // add the gene_id to the chr_to_gene[seqid] map
            } else if (rec.type == "transcript") {
                std::string gene_id = rec.attributes["Parent"];

                gene_to_transcript[rec.seqid][gene_id] = true;
                transcript_dict[rec.attributes["transcript_id"]] = Pos {rec.seqid, rec.start-1, rec.end, rec.strand[0], gene_id};
            } else if (rec.type == "exon") {
                StartEndPair sep {rec.start-1, rec.end};
                transcript_to_exon[rec.attributes["Parent"]].push_back(sep);
            }
            
            rec = parser.nextRecord();
        }

        parser.close();
    }

    // Remove duplicates from the transcript_to_exon map and also 
    // convert to Rcpp compatible data type
    List transcript_to_exon_list = remove_transcript_duplicates_to_list(transcript_to_exon, transcript_dict, false);

    
    // returning a Rcpp::List object
    List returnList = List::create(
        _["chr_to_gene"]= create_list_map_to_map(chr_to_gene), 
        _["transcript_dict"]= create_list_map_to_pos(transcript_dict), 
        _["gene_to_transcript"]= create_list_map_to_map(gene_to_transcript),
        _["transcript_to_exon"]= transcript_to_exon_list);
    std::cout << "finished parse_gff_tree\n";
    return returnList;
}

//' Parse a non gzip GTF file
List parse_gtf_tree(const char * gtf_filename) {
  // return dictionaries are: chr_to_gene, transcript_dict, gene_to_transcript, transcript_to_exon
    std::unordered_map<String, std::unordered_map<String, bool>>    chr_to_gene; // map of chr to map of genes (second map for efficient searching)
    std::unordered_map<String, Pos>                                 transcript_dict;
    std::unordered_map<String, std::unordered_map<String, bool>>    gene_to_transcript;
    std::unordered_map<String, std::list<StartEndPair>>     transcript_to_exon;

    // create the GFF3 parser
    ParseGFF3 parser (gtf_filename);

    // for every record in the GFF3 file.
    GFFRecord rec = parser.nextRecord();
    while (!parser.empty()) {
        //Rcout << "######### LINE ############\n";
        if (rec.type == "gene") {
            // add this records data to chr_to_gene
            chr_to_gene[rec.seqid][rec.attributes["gene_id"]] = true;
        } else if (rec.type == "transcript") {
            if (rec.attributes["gene_id"].length() == 0) {
                warning("Transcript did not have 'gene_id' attribute: %s", parser.formatGFFRecordAttributes(rec));
                continue;
            }

            std::string gene_id = rec.attributes["gene_id"];
            // insert the gene_id into the chr_to_gene[rec.seqid] map.
            // the chosen DS ensures there are no duplciates.
            chr_to_gene[rec.seqid][gene_id] = true;
            
            gene_to_transcript[gene_id][rec.attributes["transcript_id"]] = true;
            transcript_dict[rec.attributes["transcript_id"]] = Pos {rec.seqid, rec.start-1, rec.end, rec.strand[0], gene_id};

        } else if (rec.type == "exon") {
            if (rec.attributes["gene_id"].length() == 0) {
                warning("Exon did not have 'gene_id' attribute: %s", parser.formatGFFRecordAttributes(rec));
                continue;
            }
            chr_to_gene[rec.seqid][rec.attributes["gene_id"]] = true;
            
            // if "transcript_id" in rec.attribues keys
            if (rec.attributes["transcript_id"].length() != 0) {
                // if rec.attributes["transcript_id"] not in transcript_dict
                if (transcript_dict[rec.attributes["transcript_id"]].chr.length() == 0) {
                    gene_to_transcript[rec.attributes["gene_id"]][rec.attributes["transcript_id"]] = true;
                    transcript_dict[rec.attributes["transcript_id"]] = Pos {rec.seqid, -1, 1, rec.strand[0], rec.attributes["gene_id"]};
                }

                StartEndPair sep {rec.start-1, rec.end};
                transcript_to_exon[rec.attributes["transcript_id"]].push_back(sep);
            } else {
                warning("Exon did not have 'transcript_id' attribute: %s", parser.formatGFFRecordAttributes(rec));
            }
        }
        
         // after processing is done, get the next record
        rec = parser.nextRecord();
    }
    
    // Remove duplicates from the transcript_to_exon map and also 
    // convert to Rcpp compatible data type
    List transcript_to_exon_list = remove_transcript_duplicates_to_list(transcript_to_exon, transcript_dict, true);

    parser.close();

    // returning a Rcpp::List object
    List returnList = List::create(
        _["chr_to_gene"]= create_list_map_to_map(chr_to_gene), 
        _["transcript_dict"]= create_list_map_to_pos(transcript_dict), 
        _["gene_to_transcript"]= create_list_map_to_map(gene_to_transcript),
        _["transcript_to_exon"]= transcript_to_exon_list);
    return returnList;
}



//' Parse a GTF or GFF file
//' @description
//' Parse a Gff3 file into 3 components: chromasome to gene name, a transcript dictionary, a gene to transcript dictionary
//' and a transcript to exon dictionary.
//' These components are returned in a named list.
//'
//' @param gff_filename the file path to the gff3 file to parse
//' @return a named list with the elements 
//' "chr_to_gene", "transcript_dict", "gene_to_transcript", "transcript_to_exon", containing
//' the data parsed from the gff3 file.
//' @export
// [[Rcpp::export]]
List parse_gff_tree_cpp(const char * gff_filename) {
    // for now, we ignore gz files and gff files
    //return parse_gtf_tree(gff_filename)
    std::string filename (gff_filename);
    if (filename.find(".gtf") != std::string::npos) {
        return parse_gtf_tree(gff_filename);
    } else {
        return parse_gff_tree(gff_filename);
    }
}