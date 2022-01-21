// #include <unordered_map>
// #include <string>
// #include <fstream>
// #include <sstream>
// #include <cstring>
// #include <utility>
// #include <cctype>
// #include <Rcpp.h>

// #include "../classes/Pos.h"
// #include "ParseGFF3.hpp"
// #include "Parser.h"


// std::string ParseGFF3::formatGFFRecordAttributes(GFFRecord rec) {
//     std::stringstream stream;

//     for (std::unordered_map<std::string, std::string>::iterator it = rec.attributes.begin(); it != rec.attributes.end(); ++it) {
//         stream << it->first << "\"" << it->second << "\"" << ";";
//     }
//     return stream.str();
// }


// /// ParseGFF3 constructor. Open the file for reading
// ParseGFF3::ParseGFF3(std::string filename) {
//     // HOW CAN WE OPEN A GZ FILE FOR READING???
// 	// below is for non gzipped file
//     this->in_stream = std::ifstream (filename);
// }

// ParseGFF3::~ParseGFF3() {
// 	this->close();
// }

// bool ParseGFF3::empty(void) {
//     return this->isEmpty;
// }

// // parse a list of attributes of the form <attribute>;<attribute;... where <attribute> = <key>=<value>
// std::unordered_map<std::string, std::string> ParseGFF3::parseGTFAttributes(std::string attributes) {
// 	std::unordered_map<std::string, std::string> map {};

// 	ParseResult res = parseAttribute(attributes);
// 	ParseResult kv;
// 	while (res.first.size() != 0) {
// 		// parse the key value and add it to the map
// 		kv = parseKeyValue(res.first);
// 		if (kv.first.size() != 0) map[kv.first] = kv.second;

// 		res = parseAttribute(res.second);
// 	}

// 	// parse the final attribute
// 	kv = parseKeyValue(res.second);
// 	if (kv.first.size() != 0) map[kv.first] = kv.second;

// 	return map;
// }

// /// Parse the GFF3 file and generate a GFFRecord struct
// /// containing the required information:
// ///    seqid
// ///    source
// ///    type
// ///    start
// ///    end
// ///    score
// ///    strand
// ///    phase
// ///    attributes
// GFFRecord ParseGFF3::nextRecord() {
//     if (in_stream.is_open()) {
//         std::string line;
//         std::stringstream token;
//         std::string columns[9];

// 		getline(in_stream, line);
// 		// find the first line without a # (signifying a comment line)
// 		while (in_stream.is_open() && (line[0] == '#')) {
// 			getline(in_stream, line);
// 			continue;
// 		}

// 		// check for empty file
// 		if (line.size() != 0) {
// 			for (int i = 0; i < 8; i++) {
// 				if (line.size() == 0) {
// 					Rcpp::Rcout << "Invalid line in GFF file format.\n";
// 					return {};
// 				}

// 				ParseResult pr = parseColumn(line);
// 				columns[i] = pr.first;
// 				line = pr.second;
// 			}

// 			columns[8] = line; // fill the attributes as the final parsed end of the line

// 			return GFFRecord {
// 				(columns[0] != ".") ? columns[0] : std::string (), // seqid
// 				(columns[1] != ".") ? columns[1] : std::string (), // source
// 				(columns[2] != ".") ? columns[2] : std::string (), // type
// 				(columns[3] != ".") ? std::stoi(columns[3]) : -1, // start: int
// 				(columns[4] != ".") ? std::stoi(columns[4]) : -1, // end: int 
// 				(columns[5] != ".") ? std::stof(columns[5]) : -1, // score: float
// 				(columns[6] != ".") ? columns[6] : std::string (), // strand (+, -, .)
// 				(columns[7] != ".") ? columns[7] : std::string (), // phase
// 				(columns[8] != ".") ? this->parseGTFAttributes(columns[8]) : std::unordered_map<std::string, std::string> () // atribute unordered_map<string, string>
// 			};
// 		}
// 	}

//     /// return an empty struct if invalid input, or we reach the end of the file.
//     this->isEmpty = true; // through an error instead?
//     return {};
// }

// void ParseGFF3::close(void) {
//     in_stream.close();
// }

// std::string guess_annotation_source(std::string filename) {
//     int idx = 0;
//     // parse the file and check for "GENCODE" or "Ensembl" ???
//     // Surely there is a better way to do this
//     std::ifstream file (filename);
//     std::string line;
//     while (getline(file, line)) {
//         if (line.find("GENCODE") != std::string::npos) {
//             Rcpp::Rcout << "Parse GENCODE annotation\n";
// 			file.close();
//             return "GENCODE";
//         } else if (line.find("1\tEnsembl") != std::string::npos) {
//             Rcpp::Rcout << "Parse Ensembl annotation\n";
// 			file.close();
//             return "Ensembl";
//         }

//         if (idx++ > 1000) {
//             break;
//         }
//     }
// 	file.close();
//     return "Ensembl";
// }