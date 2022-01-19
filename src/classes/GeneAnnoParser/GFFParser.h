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
        bool            isGTF;
        std::string     annotationSource;

    public:
        GFFRecord
        parseNextRecord();
        bool
        isEmpty();
        std::string
        guessAnnotationSource();

        GFFParser(std::string filename, bool isGTF);
};

#endif