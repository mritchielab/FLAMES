#include "utility.h"

#include <string>
#include <unordered_map>
#include <vector>
#include <sstream>

#include "../classes/GeneAnnoParser/GFFRecord.h"
#include "../classes/Pos.h"

template <typename T>
void printMap(std::unordered_map<std::string, T> &m, std::function<std::string (T)> f) {
	std::stringstream ss;
	for (auto it : m) {
		ss << it.first << ": " << f(it.second) << "\n";
	}
	std::cout << ss.str();
}

std::string VecToStr(std::vector<std::string> v) {
	std::stringstream ss;
	for (auto it : v) {
		ss << it << "\t";
	}
	return ss.str();
}

std::string PosToStr(Pos p) {
	std::stringstream ss;
	ss << "Pos(chr=" << p.chr;
	ss << ", start=" << p.start;
	ss << ", end=" << p.end;
	ss << ")";
	return ss.str();
}

