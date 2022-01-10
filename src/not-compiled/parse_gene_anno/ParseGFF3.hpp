#ifndef PARSEGFF3_H
#define PARSEGFF3_H

#include <unordered_map>
#include <string>
#include <fstream>
#include <sstream>
#include <cstring>
#include <utility>
#include <cctype>
#include <Rcpp.h>
#include <iostream>

#include "Pos.h"

/// Structure for GFFRecord data used in parsing a GFF file.
/// Each struct contains the information given in a single GFF file line.
struct GFFRecord {
    std::string seqid;
    std::string source;
    std::string type;
    int start;
    int end;
    float score;
    std::string strand;
    std::string phase;
    std::unordered_map<std::string, std::string> attributes;
};

std::string guess_annotation_source(std::string);

/// Parse a GFF3 file and build a GFFRecord struct for each line.
/// GFFRecord structs are returned one at a time to provide capacity to process
/// and discard them before the rest of the file is read.
// THIS NEEDS TO BE ABLE TO READ GZ FILES
class ParseGFF3 {
    private:
        // gzFile in_stream;
        std::ifstream in_stream;
        bool isEmpty = false;
        std::unordered_map<std::string, std::string> parseGTFAttributes(std::string);
        std::unordered_map<std::string, std::string> parseGFFAttributes(std::string);
    public:
    // ignore gzip files for now
        bool empty();
        ParseGFF3 (std::string filename);
		~ParseGFF3();
        GFFRecord nextRecord(bool GFF_style_attributes=false);
        void close();
        std::string formatGFFRecordAttributes(GFFRecord rec);
};

#endif