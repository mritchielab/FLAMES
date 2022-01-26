#ifndef GENEANNOPARSER_H
#define GENEANNOPARSER_H

#include <string>
#include <unordered_map>
#include <vector>
#include <fstream>
#include <iostream>
#include <Rcpp.h>
#include <functional>

#include "GFFRecord.h"
#include "GFFParser.h"

#include "../Parser.h"
#include "../Pos.h"
#include "../StartEndPair.h"
#include "../GFFData.h"

class GeneAnnoParser;
typedef std::function<void(GFFRecord *)> ParseFunction;

class GeneAnnoParser
{
    /*
        main class to handle a GTF or GFF parser
    */

    private:
        std::string filename;
        GFFData     gffData;
        bool        isGFF;
        std::string annotationSource;

    public:
        GFFData
        parse();

        // ParseFunction
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

		static bool 
        guessGFF(std::string filename);

        GeneAnnoParser(std::string filename);
};

GFFData
parseGeneAnno(std::string filename);

bool
isFileGTF(std::string filename);

#endif // GENEANNOPARSER_H