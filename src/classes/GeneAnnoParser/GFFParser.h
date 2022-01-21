#ifndef GFFPARSER_H
#define GFFPARSER_H

#include <fstream>
#include <string>

#include "GFFRecord.h"

class GFFParser
{
    /*
        small class to read through GFF and GTF files 
    */

    private:
        std::ifstream   file;
        bool            empty;
        std::string     attributeStyle;
        std::string     annotationSource;

    public:
        GFFRecord
        parseNextRecord();
        bool
        isEmpty();

        GFFParser(std::string filename, std::string attributeStyle="GFF");
		~GFFParser();
		void close();
};

#endif