#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <filesystem>
#include <map>
#include <set>
#include <vector>
#include <algorithm>

void
gtf_to_bed(std::string in_gtf, std::string out_bed, std::string chrom_sizes_file);