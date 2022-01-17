#include <string>

#include "../classes/GFFData.h"

Rcpp::List
parse_gff_or_gtf_R(std::string filename);

GFFData
parse_gff_or_gtf(std::string filename);

GFFData
parse_gtf_tree(std::string filename);

GFFData
parse_gff_tree(std::string filename);

#endif // GFF_DATA