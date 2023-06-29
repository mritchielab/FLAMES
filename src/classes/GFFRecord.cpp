#include "GFFRecord.h"

#include <string>
#include <vector>
#include <functional>

#include "../utility/parsing.h"
#include "../utility/utility.h"

std::function<std::unordered_map<std::string, std::string>(const std::string &)> GFFRecord::chooseAttributesFunc(const std::string &filename) {
    if (filename.find(".gtf") != std::string::npos) {
        return GFFRecord::parseGTFAttributes;
    } else {
        return GFFRecord::parseGFFAttributes;
    }
}

std::unordered_map<std::string, std::string> GFFRecord::parseGFFAttributes(const std::string &attributeString) {
    if (attributeString == ".") return {};

    std::unordered_map<std::string, std::string> attributes;
    std::vector<std::string> tokens = splitStringToVector(attributeString, ';');
    for (const auto &attribute : tokens) {
        size_t commaPos = attribute.find(',');
        std::string key = attribute.substr(0, commaPos);
        std::string value = attribute.substr(commaPos);
        attributes[key] = value;
    }
    
    return attributes;
}

std::unordered_map<std::string, std::string> GFFRecord::parseGTFAttributes(const std::string &attributeString) {
    if (attributeString == ".") return {};

    std::unordered_map<std::string, std::string> attributes;
    std::vector<std::string> tokens = splitStringToVector(attributeString, ';');

    for (const auto &attribute : tokens) {
        if (attribute.size() == 0) continue;
        std::vector<std::string> items = splitStringToVector(attribute, '\"');
        if (items.size() < 2) {
            items = splitStringToVector(leftStrip(attribute), ' ');
            if (items.size() < 2) continue;
        }
        std::string key = strip(items[0]);
        std::string value = strip(items[1]);
        attributes[key] = value;
    }
    return attributes;
}



GFFRecord GFFRecord::parseGFFRecord(
    const std::string &line, 
    std::function<std::unordered_map<std::string, std::string>(const std::string &)> parseAttributesFunc) {
    
    std::vector<std::string> tokens = splitStringToVector(line, '\t');
    if (tokens.size() != 9) return {};

    GFFRecord rec;
    rec.seqname = tokens[0] == "." ? "" : tokens[0];
    rec.source = tokens[1] == "." ? "" : tokens[1];
    rec.feature = tokens[2] == "." ? "" : tokens[2];
    rec.start = stoi(tokens[3] == "." ? "-1" : tokens[3]);
    rec.end = stoi(tokens[4] == "." ? "-1" : tokens[4]);
    rec.score = stof(tokens[5] == "." ? "-1" : tokens[5]);
    rec.strand = tokens[6][0];
    rec.frame = stoi(tokens[7] == "." ? "" : tokens[7]);

    rec.attributes = parseAttributesFunc(tokens[8]);

    return rec;
}
