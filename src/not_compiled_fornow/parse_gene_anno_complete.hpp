// takes a filename
// parses it into the necessary parts

#ifndef PARSE_H
#define PARSE_H

#include <unordered_map>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <utility>


#include "Pos.h"
#include "StartEndPair.hpp"

class GFFData
{
    /*
        this is what we want to get out of the GFF file
    */

    public:
        std::unordered_map<std::string, std::vector<std::string>>
        chr_to_gene;
    
        std::unordered_map<std::string, Pos>
        transcript_dict;

        std::unordered_map<std::string, std::vector<std::string>>
        gene_to_transcript;

        std::unordered_map<std::string, std::vector<StartEndPair>>
        transcript_to_exon;

        void
        remove_transcript_duplicates(bool update_transcript_dict);
};

class Record
{
    /*
        structure for GFF or GTF record data
        each struct contains the information given in a single GFF/GTF line
    */

    private:
        std::unordered_map<std::string, std::string>
        parse_GTF_attributes(std::string attributes);

        std::unordered_map<std::string, std::string>
        parse_GFF_attributes(std::string attributes);

    public:
        std::string seqname;
        std::string source;
        std::string feature;
        int         start;
        int         end;
        float       score;
        char        strand;
        int         frame;
        std::unordered_map<std::string, std::string> attributes;

        Record(std::string line, bool is_GTF=false);
        bool
        has_attribute(std::string attribute);
        std::string
        format_attributes();
        std::string
        get_parent();
};

class Parser 
{
    /*
        a mini GTF/GFF parser for my own peace of mind
    */

    public:
        std::vector<Record>
        records;
        
        Parser(std::string filename, bool is_GTF=false);
};

GFFData
parse_GFF_or_GTF(std::string filename);

#endif