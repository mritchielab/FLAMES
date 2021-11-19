#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <filesystem>
#include <map>
#include <set>
#include <vector>
#include <algorithm>
#include <Rcpp.h>

void
gtf_to_bed_cpp(std::string in_gtf, std::string out_bed, std::string chrom_sizes_file);