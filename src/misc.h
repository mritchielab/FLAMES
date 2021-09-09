/*
  miscellaneous helper functions and things
  they will not live here forever, 
  once I figure out how to categorize them better then I will move them
*/
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>
#include <sstream>

typedef struct Pos Pos;
typedef struct Iso Iso;

std::vector<std::pair<int, int>>
pairwise (std::vector<int> input);

std::vector<int>
find_best_splice_chain(std::vector<int> raw_iso, std::vector<std::vector<int>> junction_list, int MAX_DIST);

int
if_exon_contains(std::vector<int> s1, std::vector<int> s2, int MAX_TOLERANCE);

float
get_exon_sim_pct(std::vector<int> exons1, std::vector<int> exons2);

std::vector<std::pair<int, int>>
pairwise (std::vector<int> input);

int
iv_overlap (std::pair<int, int> iv1, std::pair<int, int> iv2);

int
exon_overlap (std::vector<int> exons1, std::vector<int> exons2);

std::vector<std::pair<std::string, std::string>>
get_fa(std::string filename);
