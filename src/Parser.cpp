#include <string>
#include <fstream>
#include <functional>
#include <map>
#include <vector>
#include <utility>

#include "Parser.h"

// parse any leading spaces from the start of a token
// returns a ParseResult of {leading spaces, rest of string}
ParseResult parseSpaces(std::string full) {
	for (int i = 0; i < full.size(); i++) {
		if (!isspace(full[i])) {
			return ParseResult (full.substr(0, i), full.substr(i, std::string::npos));
		}
	}

	return ParseResult (std::string(), full);
}

// Parse a GFF3 column, each separated by spaces
// returns a ParseResult of {Column, rest of string}
ParseResult parseColumn(std::string full) {
	// parse away any leading spaces
	full = parseSpaces(full).second;
	// parse the token up until any closing spaces.
	for (int i = 0; i < full.size(); i++) {
		if (isspace(full[i])) {
			return ParseResult (full.substr(0, i), full.substr(i+1, std::string::npos));
		}
	}

	return ParseResult (full, std::string());
}
ParseResult parseColumn(std::string full, char sep) {
	// parse away any leading spaces
	full = parseSpaces(full).second;
	return parseUntilChar(full, sep);
}

std::vector<std::string> parseLine(std::string full) {
	std::vector<std::string> res;
	ParseResult p = ParseResult {std::string(), full};
	while (p.second.length()) {
	    p = parseColumn(p.second);
	    res.push_back(p.first);
	}
	return res;
}
std::vector<std::string> parseLine(std::string full, char sep) {
	std::vector<std::string> res;
	ParseResult p = ParseResult {std::string(), full};
	while (p.second.length()) {
		p = parseColumn(p.second, sep);
		res.push_back(p.first);
	}
	return res;
}

std::map<std::string, int> parsePairsToMap(std::ifstream &input) {
	std::map<std::string, int> map;
	std::string line;
	while (std::getline(input, line)) {
		std::vector<std::string> tokens = parseLine(line);
		map[tokens[0]] = std::stoi(tokens[1]);
	}

	return map;
}

// parse a key value pair separated by =
// returns a ParseResult of {key, value}
ParseResult parseKeyValue(std::string full) {
	return parseKeyValue(full, '=');
}
ParseResult parseKeyValue(std::string full, char seperator) {
	int sep = full.find(seperator);
	if (sep != std::string::npos) {
		return ParseResult (full.substr(0, sep), full.substr(sep+1, std::string::npos));
	}

	return ParseResult ();
}

// parse an attribute from a colon separated list
// returns a ParseResult of {attribute, rest of string}
ParseResult parseAttribute(std::string full) {
	// parse away any leading spaces
	full = parseSpaces(full).second;
	int colonPos = full.find(';');

	if (colonPos != std::string::npos) {
		return ParseResult (full.substr(0, colonPos), full.substr(colonPos+1, std::string::npos));
	}

	return ParseResult (full, std::string());
}