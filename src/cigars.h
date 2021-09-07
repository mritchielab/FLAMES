#include <string>
#include <vector>

std::string
generate_cigar (std::vector <std::pair <int, int>> cigar);

std::vector<std::pair<int, int>>
smooth_cigar (std::vector<std::pair<int, int>> cigar, int threshold=10);