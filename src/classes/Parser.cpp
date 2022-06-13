#include "Parser.h"

#include <string>
#include <fstream>
#include <functional>
#include <map>
#include <vector>
#include <utility>
#include <cctype>

// parse any leading spaces from the start of a token
// returns a ParseResult of {leading spaces, rest of string}
ParseResult parseSpaces(std::string full) {
	return parseLeadingChar(full, isspace);
}

// parse any leading char away from the start of a token
// returns a ParseResult of {leading chars, rest of string}
ParseResult parseLeadingChar(std::string full, char c) {
	for (int i = 0; i < (int)full.size(); ++i) {
		// if we've encountered a character not of the given char, return the remaining strin
		if (full[i] != c) {
			return ParseResult (full.substr(0, i), full.substr(i, std::string::npos));
		}
	}

	return ParseResult (full, std::string());
}
ParseResult parseLeadingChar(std::string full, std::function<int(int)> condition) {
	for (int i = 0; i < (int)full.size(); ++i) {
		if (!condition(full[i])) {
			return ParseResult(full.substr(0, i), full.substr(i, std::string::npos));
		}
	}

	return ParseResult(full, std::string());
}

// Parse a GFF3 column, each separated by spaces
// returns a ParseResult of {Column, rest of string}
ParseResult parseColumn(std::string full) {
	// parse away any leading spaces
	full = parseSpaces(full).second;
	// parse the token up until any closing spaces.
	for (int i = 0; i < (int)full.size(); ++i) {
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
	if (sep != (int)std::string::npos) {
		return ParseResult (full.substr(0, sep), full.substr(sep+1, std::string::npos));
	}

	return ParseResult ();
}

// GTF attributes follow a specific convention:
// key "val"; key "val"; another key "different value"
ParseResult parseGTFKeyValue(std::string full) {
	// parse away leading spaces
	ParseResult res = parseSpaces(full);
	// just get the key, it should be followed by a space
	res = parseUntilChar(res.second, ' ');
	auto key = res.first;
	// parse away next space
	res = parseSpaces(res.second);
	// then parse away the first quote
	res = parseLeadingChar(res.second, '\"');
	// then parse until the second quote
	res = parseUntilChar(res.second, '\"');
	auto val = res.first;
	return {key, val};
}

// parse an attribute from a semicolon separated list
// returns a ParseResult of {attribute, rest of string}
ParseResult parseAttribute(std::string full) {
	// parse away any leading spaces
	full = parseSpaces(full).second;

	int colonPos = full.find(';');

	if (colonPos != (int)std::string::npos) {
		return ParseResult (full.substr(0, colonPos), full.substr(colonPos+1, std::string::npos));
	}

	return ParseResult (full, std::string());
}