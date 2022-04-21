#include "utility.h"

#include <string>
#include <unordered_map>
#include <vector>
#include <sstream>

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

StartEndPair mostCommonSEP(const std::vector<StartEndPair> &values) {
	StartEndPair newValue = values[0];
	std::unordered_map<int, int> startValues, endValues;
	for (const auto &it : values) {
		startValues[it.start]++;
		endValues[it.end]++;
	}

	int startMaxCount = 0;
	int endMaxCount = 0;
	for (const auto &[start, count] : startValues) {
		if (count > startMaxCount) {
			startMaxCount = count;
			newValue.start = start;
		} else if (count == startMaxCount && start < newValue.start) {
			newValue.start = start;
		}
	}
	for (const auto &[end, count] : endValues) {
		if (count > endMaxCount) {
			endMaxCount = count;
			newValue.end = end;
		} else if (count == endMaxCount && end > newValue.end) {
			newValue.end = end;
		}
	}

	return newValue;
}
