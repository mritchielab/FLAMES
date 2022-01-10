#ifndef PARSER_H
#define PARSER_H

#include <string>
#include <functional>
#include <map>
#include <utility>
#include <vector>
#include <iostream>
#include <unordered_map>

typedef std::pair<std::string, std::string> ParseResult;

// parse any leading spaces from the start of a token
// returns a ParseResult of {leading spaces, rest of string}
ParseResult parseSpaces(std::string);
ParseResult parseLeadingChar(std::string, char);

// Parse a GFF3 column, each separated by spaces
// returns a ParseResult of {Column, rest of string}
ParseResult parseColumn(std::string);
ParseResult parseColumn(std::string, char);

std::vector<std::string> parseLine(std::string);
std::vector<std::string> parseLine(std::string, char);

inline ParseResult parseUntilChar(std::string full, char tok) {
	int tok_pos = full.find(tok);
	if (tok_pos != std::string::npos) {
		return ParseResult {full.substr(0, tok_pos), full.substr(tok_pos + 1, full.length())};
	} else {
		return ParseResult {full, std::string()};
	}
}

std::map<std::string, int> parsePairsToMap(std::ifstream &);

// parse a key value pair separated by =
// returns a ParseResult of {key, value}
ParseResult parseKeyValue(std::string);
ParseResult parseKeyValue(std::string, char);
ParseResult parseGTFKeyValue(std::string);

// parse an attribute from a colon separated list
// returns a ParseResult of {attribute, rest of string}
ParseResult parseAttribute(std::string);
std::unordered_map<std::string, std::string>
parseGTFAttributes(std::string);
#endif // PARSER_H