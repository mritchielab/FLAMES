#ifndef GFFRECORD_H
#define GFFRECORD_H

#include <unordered_map>
#include <string>
#include <vector>
#include <sstream>

#include "../Parser.h"

typedef std::unordered_map<std::string, std::string> AttributesMap;

class GFFRecord
{
    /*
        structure for GFF or GTF record data
        each struct contains the information given in a single GFF/GTF line
    */

    private:
        AttributesMap
        parseAttributes(std::string attributes, std::string attributeStyle);

        AttributesMap
        parseGTFAttributes(std::string attributes);

        AttributesMap
        parseGFFAttributes(std::string attributes);

    public:
        std::string     seqname;
        std::string     source;
        std::string     feature;
        int             start;
        int             end;
        float           score;
        char            strand;
        int             frame;
        AttributesMap   attributes;

        bool            broken; // for when the record is just a placeholder for an illegal line

        bool
        hasAttribute(std::string attribute);

        std::string
        printAttributes();

        std::string
        print();

        GFFRecord(std::string line, std::string attributeStyle="GFF");
        GFFRecord();
};

#endif