#ifndef GFFRECORD_H
#define GFFRECORD_H

#include <unordered_map>
#include <string>
#include <vector>
#include <functional>

struct GFFRecord {
    std::string seqname, source, feature;
    int start, end;
    float score;
    char strand;
    int frame;
    std::unordered_map<std::string, std::string> attributes;

    GFFRecord() {};
    
    static GFFRecord parseGFFRecord(
        const std::string &line, 
        std::function<std::unordered_map<std::string, std::string>(const std::string &)> parseAttributesFunc);

    // static std::function<GFFRecord(const std::string &)> parseGFFRecordFactoryFromAttributes(const std::string &filename);
    static std::unordered_map<std::string, std::string> parseGFFAttributes(const std::string &attributes);
    static std::unordered_map<std::string, std::string> parseGTFAttributes(const std::string &attributes);
    static std::function<std::unordered_map<std::string, std::string>(const std::string &)> chooseAttributesFunc(const std::string &filename);

};

#endif // GFFRECORD_H