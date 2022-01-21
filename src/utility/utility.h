#include <string>
#include <unordered_map>
#include <vector>

#include "../classes/GeneAnnoParser/GFFRecord.h"
#include "../classes/Pos.h"

template <typename T>
void printMap(std::unordered_map<std::string, T> &, std::function<std::string (T)>);

std::string PosToStr(Pos);

std::string VecToStr(std::vector<std::string>);