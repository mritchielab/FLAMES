#include "GFFRecord.h"

#include <unordered_map>
#include <string>
#include <vector>
#include <sstream>

#include "../Parser.h"

/*
    takes a line of GFF or GTF,
    parses it into Record format
*/
GFFRecord::GFFRecord
(
    std::string line, 
    std::string attributeStyle
)
{
    // if the line is a comment, the record is considered broken
    if (line.size() >= 1 && line[0] == '#') {
        this->broken = true;
        return;
    }
     
    auto
    columns = parseLine(line, '\t');

    // if the line is not spaced correctly, the record is considered broken
    if (columns.size() != 9) {
        this->broken = true;
        return;
    }
    
    this->seqname   = (columns[0] != ".") ? columns[0] : "";
    this->source    = (columns[1] != ".") ? columns[1] : "";
    this->feature   = (columns[2] != ".") ? columns[2] : "";
    this->start     = (columns[3] != ".") ? std::stoi(columns[3]) : -1;
    this->end       = (columns[4] != ".") ? std::stoi(columns[4]) : -1;
    this->score     = (columns[5] != ".") ? std::stof(columns[5]) : -1;
    this->strand    = (columns[6] != ".") ? columns[6][0] : '+';
    this->frame     = (columns[7] != ".") ? std::stoi(columns[7]) : -1;
    this->attributes= (columns[8] != ".") ? GFFRecord::parseAttributes(columns[8], attributeStyle) : AttributesMap();

    this->broken = false;
}


/*
    initialises a broken record
*/
GFFRecord::GFFRecord()
{
    this->broken = true;
}

/*
    takes the attributes section of a GFF or GTF line,
    selects the right parsing method and parses it
*/
AttributesMap
GFFRecord::parseAttributes(std::string attributes, std::string attributeStyle)
{
    return attributeStyle == "GFF" ? 
        parseGFFAttributes(attributes) : parseGTFAttributes(attributes);
}

/*  
    GTF files follow a different attribute convention
    this parser handles them
    (the second format on this site: https://m.ensembl.org/info/website/upload/gff.html)
*/
AttributesMap
GFFRecord::parseGTFAttributes(std::string attributes)
{
    AttributesMap attributesMap;

	ParseResult res = parseAttribute(attributes);
	ParseResult keyValue;
	while (res.first.size() != 0) {
		// parse the key value and add it to the map
        ParseResult keyValue = parseGTFKeyValue(res.first);
		if (keyValue.first.size() != 0) attributesMap[keyValue.first] = keyValue.second;

		res = parseAttribute(res.second);
	}

	// parse the final attribute
	keyValue = parseGTFKeyValue(res.second);
	if (keyValue.first.size() != 0) attributesMap[keyValue.first] = keyValue.second;
    
    return attributesMap;
}

/*
    takes the attributes in GFF3 format,
    parses them into an AttributesMap by splitting on ';' and '='
*/
AttributesMap
GFFRecord::parseGFFAttributes(std::string attributes)
{
    AttributesMap attributesMap;

    // separate out each attribute
    auto attributesVector = parseLine(attributes, ';');
    for (const auto & attribute : attributesVector) {
        // then parse each one into key and value
        auto keyValue = parseKeyValue(attribute);
        // and put them in the map
        attributesMap[keyValue.first] = keyValue.second;
    }

    return attributesMap;
}

/*
    checks whether this record has a particular attribute
*/
bool
GFFRecord::hasAttribute(std::string attribute)
{
    return this->attributes.count(attribute) > 0;
}

/*
    prints out the attributes of a record to string format
*/
std::string
GFFRecord::printAttributes()
{
	std::stringstream ss;
	for (auto it : this->attributes) {
		ss << it.first << "=" << it.second << ";";
	}
	return ss.str();
}

/*
    prints out the entire record to a string format
*/
std::string
GFFRecord::print()
{
	std::string attr = GFFRecord::printAttributes();

	std::stringstream ss;
	ss << "GFFRecord(seqname=" << this->seqname
		<< ", source=" << this->source
		<< ", feature=" << this->feature
		<< ", start=" << this->start
		<< ", end=" << this->end 
		<< ", score=" << this->score 
		<< ", strand=" << this->strand
		<< ", frame=" << this->frame
		<< ", attributes={" << attr << "}";
	return ss.str();
}