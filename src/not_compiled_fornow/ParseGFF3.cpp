
#include "Pos.h"
#include "ParseGFF3.hpp"
#include "Parser.hpp"


std::string ParseGFF3::formatGFFRecordAttributes(GFFRecord rec) {
    std::stringstream stream;

    for (std::unordered_map<std::string, std::string>::iterator it = rec.attributes.begin(); it != rec.attributes.end(); ++it) {
        stream << it->first << "\"" << it->second << "\"" << ";";
    }
    return stream.str();
}


/// ParseGFF3 constructor. Open the file for reading
ParseGFF3::ParseGFF3(std::string filename) {
    // HOW CAN WE OPEN A GZ FILE FOR READING???
	// below is for non gzipped file
    std::cout << "created Parser\n";
    this->in_stream = std::ifstream (filename);
}

ParseGFF3::~ParseGFF3() {
	this->close();
}

bool ParseGFF3::empty(void) {
    return this->isEmpty;
}

// parse a list of attributes of the form <attribute>;<attribute;... where <attribute> = <key>=<value>
std::unordered_map<std::string, std::string>
ParseGFF3::parseGFFAttributes(std::string attributes)
{
    std::unordered_map<std::string, std::string> map {};

	ParseResult res = parseAttribute(attributes);
	ParseResult kv;
	while (res.first.size() != 0) {
		// parse the key value and add it to the map
		kv = parseKeyValue(res.first);
		if (kv.first.size() != 0) map[kv.first] = kv.second;

		res = parseAttribute(res.second);
	}

	// parse the final attribute
	kv = parseKeyValue(res.second);
	if (kv.first.size() != 0) map[kv.first] = kv.second;

    for (const auto & [k, v] : map) {
        std::cout << "\t" << k << ":" << v << "\n";
    }
    return map;
}


std::unordered_map<std::string, std::string>
ParseGFF3::parseGTFAttributes
(
    std::string attributes
)
{
    /*  
        GTF files follow a different attribute convention
        this parser handles them
        (the second format on this site: https://m.ensembl.org/info/website/upload/gff.html)
    */
   	std::cout << "started parseGTFAttributes\n";
    std::unordered_map<std::string, std::string> map {};

	ParseResult res = parseAttribute(attributes);
	ParseResult kv;
	while (res.first.size() != 0) {
		// parse the key value and add it to the map
        ParseResult kv = parseGTFKeyValue(res.first);
		if (kv.first.size() != 0) map[kv.first] = kv.second;

		res = parseAttribute(res.second);
	}

	// parse the final attribute
	kv = parseGTFKeyValue(res.second);
	if (kv.first.size() != 0) map[kv.first] = kv.second;
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
GFFRecord ParseGFF3::nextRecord(bool GFF_style_attributes) {
    std::cout << "started nextRecord\n";
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
                std::cout << "Invalid line in GFF file format.\n";
                return {};
            }

            
            // extract the attributes, using the appropriate style
            auto
            attributes_map = GFF_style_attributes ?
                this->parseGFFAttributes(columns[8]) : this->parseGTFAttributes(columns[8]);
            std::cout << "created record from line:" << line << "\n";
            return GFFRecord {
                (columns[0] != ".") ? columns[0] : std::string (), // seqid
                (columns[1] != ".") ? columns[1] : std::string (), // source
                (columns[2] != ".") ? columns[2] : std::string (), // type
                (columns[3] != ".") ? std::stoi(columns[3]) : -1, // start: int
                (columns[4] != ".") ? std::stoi(columns[4]) : -1, // end: int 
                (columns[5] != ".") ? std::stof(columns[5]) : -1, // score: float
                (columns[6] != ".") ? columns[6] : std::string (), // strand (+, -, .)
                (columns[7] != ".") ? columns[7] : std::string (), // phase
                (columns[8] != ".") ? attributes_map : std::unordered_map<std::string, std::string> () // attribute unordered_map<string, string>
            };
        }
        this->isEmpty = true; // through an error instead?
	}

    /// return an empty struct if invalid input, or we reach the end of the file.
    this->isEmpty = true; // through an error instead?
    return {};
}

void ParseGFF3::close(void) {
    in_stream.close();
}

std::string guess_annotation_source(std::string filename) {
    int idx = 0;
    // parse the file and check for "GENCODE" or "Ensembl" ???
    // Surely there is a better way to do this
    std::ifstream file (filename);
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
