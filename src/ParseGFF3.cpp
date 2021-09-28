#include "ParseGFF3.hpp"

std::string ParseGFF3::formatGFFRecordAttributes(GFFRecord rec) {
    std::stringstream stream;

    for (std::unordered_map<std::string, std::string>::iterator it = rec.attributes.begin(); it != rec.attributes.end(); ++it) {
        stream << it->first << "\"" << it->second << "\"" << ";";
    }
    return stream.str();
}


/// ParseGFF3 constructor. Open the file for reading
ParseGFF3::ParseGFF3(const char * filename) {
    // HOW CAN WE OPEN A GZ FILE FOR READING???
    this->in_stream = std::ifstream (filename);
}

bool ParseGFF3::empty(void) {
    return this->isEmpty;
}

std::unordered_map<std::string, std::string> ParseGFF3::parseGTFAttributes(std::string attributes) {
    std::unordered_map<std::string, std::string> map {};

    std::stringstream token;
    std::string key;
    bool firstQ = true;
    for (std::string::iterator it = attributes.begin(); it != attributes.end(); ++it) {
        if (*it == ';') {
            map[key] = token.str();
            token.str(std::string ());
            firstQ = true;
        } else if (*it == '\"') {
            if (firstQ) {
                key = token.str();
                token.str(std::string());
                firstQ = false;
            }
        } else if (*it != ' ') {
            token << *it;
        }
    }
    // insert the final pair into the map
    //map[key] = token.str();

    return map;
}

/// Parse the GFF3 file and generate a GFFRecord struct
/// containing the required information:
///    seqid
///    source
///    type
///    start
///    end
///    score
///    strand
///    phase
///    attributes
GFFRecord ParseGFF3::nextRecord() {
    if (in_stream.is_open()) {
        std::string line;
        std::stringstream token;
        std::string columns[9];
        int cur_col = 0;

        if (getline(in_stream, line)) {
            // iterate over the current and grab the 9 tab separated columns
            for (std::string::iterator it = line.begin(); it != line.end(); ++it) {
                // loop through every character in the line
                if (*it == '\t') {
                    // insert the token into the column array and increment column
                    columns[cur_col++] = token.str();
                    token.str(std::string ()); // reset the stringstream by setting it to the empty string
                } else if (*it != ' ') {
                    // insert the character into the current token if not whitespace
                    token << *it;
                }
            }
            // after end of line need to flush the stringstream to the array
            columns[cur_col++] = token.str();

            // make sure GFF line is valid
            if (cur_col != 9) {
                Rcout << "Invalid line in GFF file format.\n";
                return {};
            }

            return GFFRecord {
                (columns[0] != ".") ? columns[0] : std::string (), // seqid
                (columns[1] != ".") ? columns[1] : std::string (), // source
                (columns[2] != ".") ? columns[2] : std::string (), // type
                (columns[3] != ".") ? std::stoi(columns[3]) : -1, // start: int
                (columns[4] != ".") ? std::stoi(columns[4]) : -1, // end: int 
                (columns[5] != ".") ? std::stof(columns[5]) : -1, // score: float
                (columns[6] != ".") ? columns[6] : std::string (), // strand (+, -, .)
                (columns[7] != ".") ? columns[7] : std::string (), // phase
                (columns[8] != ".") ? this->parseGTFAttributes(columns[8]) : std::unordered_map<std::string, std::string> () // atribute unordered_map<string, string>
            };
        }

        this->isEmpty = true;
        return {};
    }

    /// return an empty struct if invalid input, or we reach the end of the file.
    this->isEmpty = true; // through an error instead?
    return {};
}

void ParseGFF3::close(void) {
    in_stream.close();
}

bool StartEndPairCompare(const StartEndPair &a, const StartEndPair &b) {
    // compare a and b, return true if a is 'less than' b
    // in this case, 'less than' is defined if a.start is less than b.start
    return a.start < b.start;
}

std::string guess_annotation_source(const char * filename) {
    int idx = 0;
    // parse the file and check for "GENCODE" or "Ensembl" ???
    // Surely there is a better way to do this
    std::ifstream file (filename);
    std::string line;
    while (getline(file, line)) {
        if (line.find("GENCODE") != std::string::npos) {
            Rcout << "Parse GENCODE annotation\n";
            return "GENCODE";
        } else if (line.find("1\tEnsembl") != std::string::npos) {
            Rcout << "Parse Ensembl annotation\n";
            return "Ensembl";
        }

        if (idx++ > 1000) {
            break;
        }
    }

    return "Ensembl";
}