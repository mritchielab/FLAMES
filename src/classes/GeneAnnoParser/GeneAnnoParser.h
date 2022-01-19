#ifndef GENEANNOPARSER_H
#define GENEANNOPARSER_H

#include <string>
#include <unordered_map>
#include <vector>
#include <fstream>
#include <iostream>
#include <Rcpp.h>

#include "GFFRecord.h"
#include "GFFParser.h"

#include "../Parser.h"
#include "../Pos.h"
#include "../StartEndPair.h"
#include "../GFFData.h"

class GeneAnnoParser;
typedef void (GeneAnnoParser::*ParseFunction)(GFFRecord*);

class GeneAnnoParser
{
    /*
        main class to handle a GTF or GFF parser
    */

    private:
        std::string filename;
        GFFData     gffData;
        bool        isGTF;
        std::string annotationSource;
        GFFParser * gffParser;

    public:
        GFFData
        parse();

        ParseFunction
        selectParseFunction();

        void
        parseGTF(GFFRecord * rec);

        void
        parseGFF(GFFRecord * rec);

        void
        parseEnsembl(GFFRecord * rec);

        void
        parseGENCODE(GFFRecord * rec);

        GeneAnnoParser(std::string filename);
        ~GeneAnnoParser();
};

GFFData
parseGeneAnno(std::string filename);

bool
isFileGTF(std::string filename);

#endif // GENEANNOPARSER_H