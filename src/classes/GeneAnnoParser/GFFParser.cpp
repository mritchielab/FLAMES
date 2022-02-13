#include "GFFParser.h"

#include <fstream>
#include <string>

#include <Rcpp.h>

#include "GFFRecord.h"

/*
    initialises the GFFParser, opening the file
	attributeStyle can be either GFF or GTF (this should really be a bool flag)
*/
GFFParser::GFFParser(std::string filename, bool isGFF)
{
    // open the file
    this->file  = std::ifstream(filename);
    this->empty = false;
    this->isGFF = isGFF; 
}

/*
    parses the next record in the file,
    returning a completed GFFRecord
*/
GFFRecord
GFFParser::parseNextRecord()
{
    if (file.is_open()) {
        std::string line;
        if (getline(file, line)) {
			while (line[0] == '#')	 {
				getline(file, line);
			}
            // turn each line into a record
            GFFRecord rec(line, isGFF);
            return rec;
        } else {
			// if we are out of lines, return nothing
			this->empty = true;
			return GFFRecord();
		}
    }
    Rcpp::Rcout << "file is not open\n";
    return GFFRecord();
}

/*
    checks whether the parser is empty
*/
bool
GFFParser::isEmpty()
{
    return this->empty;
}

GFFParser::~GFFParser() {
	this->close();
}

void GFFParser::close() {
	file.close();
}
/*
    parses the file and checks for "GENCODE" or "Ensembl"
    surely there's a better way to do this
*/
std::string 
GFFParser::guessAnnotationSource()
{
    int idx = 0;
    std::string line;
    while (getline(file, line)) {
        if (line.find("GENCODE") != std::string::npos) {
			file.close();
            return "GENCODE";
        } else if (line.find("1\tEnsembl") != std::string::npos) {
			file.close();
            return "Ensembl";
        }

        if (idx++ > 1000) {
            break;
        }
    }
	file.clear();
	file.seekg(0);
    return "Ensembl";
}