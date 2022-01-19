#include "GFFParser.h"

/*
    initialises the GFFParser, opening the file
*/
GFFParser::GFFParser(std::string filename, bool isGTF)
{
    // open the file
    this->file  = std::ifstream(filename);
    this->empty = false;
    this->isGTF = isGTF;
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
            std::cout << line << "\n";
            if (line[0] == '#') {
                return GFFRecord();
            }
            // turn each line into a record
            GFFRecord * rec = new GFFRecord(line, this->isGTF);
            return *rec;
        }
        // if we are out of lines, return nothing
        this->empty = true;
        return GFFRecord();
    }
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
	file.close();
    return "Ensembl";
}