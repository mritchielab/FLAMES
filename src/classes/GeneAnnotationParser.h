#ifndef GENE_ANNOTATION_PARSER_H
#define GENE_ANNOTATION_PARSER_H

#include <string>

#include "GFFData.h"

GFFData parse_gff_file(const std::string &gff_filename);

#endif // GENE_ANNOTATION_PARSER_H