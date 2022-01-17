#include "GFFParser.h"

/*
    initialises the GFFParser, opening the file
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
    std::cout << "started parseNextRecord\n";
    if (file.is_open()) {
        std::cout << "file is open\n";
        std::string line;
        if (getline(file, line)) {
            std::cout << "gotline successfully, and it was:\n";
            std::cout << line << "\n";
            if (line[0] == '#') {
                std::cout << "skipping commment\n";
                return GFFRecord();
            }
            // turn each line into a record
            GFFRecord * rec = new GFFRecord(line, this->attributeStyle);
            return *rec;
        }
        // if we are out of lines, return nothing
        this->empty = true;
        std::cout << "couldn't getline, file isEmpty\n";
        return GFFRecord();
    }
    std::cout << "file is not open\n";
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
            std::cout << "Parse GENCODE annotation\n";
			file.close();
            return "GENCODE";
        } else if (line.find("1\tEnsembl") != std::string::npos) {
            std::cout << "Parse Ensembl annotation\n";
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