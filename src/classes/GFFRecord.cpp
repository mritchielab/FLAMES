#include "GFFRecord.h"

#include <string>
#include <vector>
#include <functional>

#include "../utility/parsing.h"
#include "../utility/utility.h"

// Create a function which will parse a GFFRecord using a chosen attributes parsing function
// Return that function to the caller, capturing the attributes function
// std::function<GFFRecord(const std::string &)> GFFRecord::parseGFFRecordFactoryFromAttributes(const std::string &filename) {
//     auto attributesFunc = 
//         GFFRecord::chooseAttributesFunc(filename);
//     return [&attributesFunc](const std::string &line) {
//         return GFFRecord::parseGFFRecord(line, attributesFunc);
//     };
// }

std::function<std::unordered_map<std::string, std::string>(const std::string &)> GFFRecord::chooseAttributesFunc(const std::string &filename) {
    if (filename.find(".gtf") != std::string::npos) {
        return GFFRecord::parseGTFAttributes;
    } else {
        return GFFRecord::parseGFFAttributes;
    }
}

std::unordered_map<std::string, std::string> GFFRecord::parseGFFAttributes(const std::string &attributeString) {
    if (attributeString == ".") return {};

    std::vector<std::string> tokens = splitStringToVector(strip(attributeString), ';');
    std::unordered_map<std::string, std::string> attributes;
    
    for (const std::string &attribute : tokens) {
        size_t commaPos = attribute.find('=');
        std::string key = attribute.substr(0, commaPos);
        std::string value = attribute.substr(commaPos + 1);
        attributes[key] = value;
    }

    return attributes;
}

std::unordered_map<std::string, std::string> GFFRecord::parseGTFAttributes(const std::string &attributeString) {
    if (attributeString == ".") return {};

    std::vector<std::string> tokens = splitStringToVector(attributeString, ';');
    std::unordered_map<std::string, std::string> attributes;

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



GFFRecord GFFRecord::parseGFFRecord(const std::string &line, 
        std::function<std::unordered_map<std::string, std::string>(const std::string &)> parseAttributesFunc) {
    
    std::vector<std::string> tokens = splitStringToVector(line, '\t');
    if (tokens.size() != 9) return {};

    GFFRecord rec;
    rec.seqname = tokens[0] == "." ? "" : tokens[0];
    rec.source = tokens[1] == "." ? "" : tokens[1];
    rec.feature = tokens[2] == "." ? "" : tokens[2];
    rec.start = tokens[3] != "." ? parseInt(tokens[3]) : -1;
    rec.end = tokens[4] != "." ? parseInt(tokens[4]) : -1;
    rec.score = tokens[5] != "." ? parseFloat(tokens[5]) : -1.0;
    rec.strand = tokens[6][0];
    rec.frame = tokens[7] != "." ? parseInt(tokens[7]) : -1;

    rec.attributes = parseAttributesFunc(tokens[8]);
    return rec;
}
