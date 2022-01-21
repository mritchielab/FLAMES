#include "GFFParser.h"

#include <fstream>
#include <string>

#include <Rcpp.h>

#include "GFFRecord.h"

/*
    initialises the GFFParser, opening the file
	attributeStyle can be either GFF or GTF (this should really be a bool flag)
*/
GFFParser::GFFParser(std::string filename, std::string attributeStyle)
{
    // open the file
    this->file = std::ifstream(filename);
    this->empty = false;
    this->attributeStyle = attributeStyle;
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
            GFFRecord rec(line, this->attributeStyle);
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