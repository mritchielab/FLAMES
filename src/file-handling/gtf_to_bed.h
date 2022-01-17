#ifndef GTF_TO_BED_H
#define GTF_TO_BED_H

#include <string>
#include <sstream>
#include <fstream>
#include <filesystem>
#include <map>
#include <set>
#include <vector>
#include <algorithm>
#include <functional>
#include <numeric>
#include <Rcpp.h>

#include "../classes/Parser.h"

void
gtf_to_bed(std::string in_gtf, std::string out_bed, std::string chrom_sizes_file);

#endif // GTF_TO_BED_H