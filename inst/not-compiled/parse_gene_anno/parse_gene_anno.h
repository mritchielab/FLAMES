#include <unordered_map>
#include <forward_list>
#include <list>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <Rcpp.h>

#include "ParseGFF3.hpp"
#include "StartEndPair.hpp"

using namespace Rcpp;

List 
create_list_map_to_map(std::unordered_map<String, std::unordered_map<String, bool>> &map);

List
create_list_map_to_pos(std::unordered_map<String, Pos> &map);

List
remove_transcript_duplicates_to_list(std::unordered_map<String, std::list<StartEndPair>> &transcript_to_exon, std::unordered_map<String, Pos> &transcript_dict, bool update_transcript_dict);

List
parse_gff_tree(const char * gff_filename);

List
parse_gtf_tree(const char * gtf_filename);

List
parse_gff_tree_cpp(const char * gff_filename);