#ifndef UTILITY_H
#define UTILITY_H

#include <string>
#include <sstream>
#include <unordered_map>
#include <vector>
#include <functional>
#include <algorithm>
#include <set>
#include <utility>
#include <cctype>

#include <Rcpp.h>

#include "../classes/Pos.h"
#include "../classes/StartEndPair.h"

namespace ranges
{
    template<typename Range, typename Function>
    Function for_each(Range& range, Function f)
    {
        return std::for_each(std::begin(range), std::end(range), f);
    }

	template<typename Range, typename Function>
	Function transform(Range& range, Range& out, Function f) {
		return std::transform(std::begin(range), std::end(range), std::begin(out), f);
	}

	template<typename T, typename U>
	std::vector<U> map(const std::vector<T> &range, std::function<U(const T &)> f) {
		std::vector<U> out;
		for (const T &t : range) {
			out.push_back(f(t));
		}
		return out;
	}

	template <typename T, typename U>
	U sumMap(const std::vector<T> &vec, std::function<U(T)> f) {
		U sum = 0;
		for (auto it : vec) {
			sum += f(it);
		}
		return sum;
	}

	template <typename T>
	void filter(const std::vector<T> &vec, std::vector<T> &out, std::function<bool(T)> f) {
		std::copy_if(vec.begin(), vec.end(), std::back_inserter(out), f);
	}

	template <typename T>
	std::vector<T> filter(const std::vector<T> &vec, std::function<bool(T)> f) {
		std::vector<T> out;
		filter(vec, out, f);
		return out;
	}
	
	template <typename T>
	T slice(const T &vec, int start, int end) {
		return {std::next(vec.begin(), start), std::next(vec.begin(), end)};
	}

	template <typename T, typename U>
	int count(const T &vec, const U &val) {
		return std::count(vec.begin(), vec.end(), val);
	}

	template <typename T>
	bool hasDuplicates(const std::vector<T> &vec) {
		std::set<T> set(vec.begin(), vec.end());
		return vec.size() > set.size();
	}

	template <typename T>
	void sort(std::vector<T> &vec) {
		std::sort(vec.begin(), vec.end());
	}
	template <typename T, class Comp>
	void sort(std::vector<T> &vec, Comp f) {
		std::sort(vec.begin(), vec.end(), f);
	}
	
	template <typename T>
	bool doesNotContain(const std::vector<T> &vec, const T &val) {
		return std::find(vec.begin(), vec.end(), val) == vec.end();
	}
}


template <typename T>
void printMap(std::unordered_map<std::string, T> &m, std::function<std::string (T)> f) {
	std::stringstream ss;
	for (auto it : m) {
		ss << it.first << ": " << f(it.second) << "\n";
	}
	std::cout << ss.str();
}

inline std::string PosToStr(const Pos &p) {
	std::stringstream ss;
	ss << "Pos(chr=" << p.chr;
	ss << ", start=" << p.start;
	ss << ", end=" << p.end;
	ss << ")";
	return ss.str();
}

// inline std::string VecToStr(const std::vector<std::string> &v) {
// 	std::stringstream ss;
// 	for (auto it : v) {
// 		ss << it << "\t";
// 	}
// 	return ss.str();
// }

template <typename T>
inline std::string VecToStr(const std::vector<T> &v, char open='(', char close=')') {
	std::stringstream ss;
	ss << open;
	if (v.size() > 0) {
		ss << v.at(0);
		for (size_t i = 1; i < v.size(); i++) {
			ss << ", " << v.at(i);
		}
	}
	ss << close;
	return ss.str();
}
/*
 * At vector position, determine the most value across all vectors
 * @return a vector of the most common elements at each position
*/
template <typename T>
inline std::vector<T> 
mostCommonEachCell(const std::vector<std::vector<T>> &values, int maxSize) {
	std::vector<T> newValue;
	for (int i = 0; i < maxSize; i++) {
		// thisPosValueCounts is a map of every value in the column we're inspecting
		// with it's corresponding appearance count.
		std::unordered_map<T, int> thisPosValueCounts;
		for (int j = 0, size = values.size(); j < size; j++) {
			thisPosValueCounts[values[j][i]]++;
		}
		// get the the most common value of junction
		// at each i position 
		T bestValue;
		int maxCount = 0;
		for (auto &[value, count] : thisPosValueCounts) {
			if (count > maxCount) {
				bestValue = value;
				maxCount = count;
			}
		}
		newValue.push_back(bestValue);
	}

	return newValue;
}

template <typename T>
inline T mostCommon(const std::vector<T> &values) {
	std::unordered_map<T, int> counts;
	for (auto it : values) {
		counts[it]++;
	}

	T max = values[0];
	for (const auto &[val, count] : counts) {
		max = count > counts[max] ? val : max;
	}

	return max;
}

inline StartEndPair mostCommonSEP(const std::vector<StartEndPair> &values) {
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

/*
 * Taking a vector of values (which contains duplicates)
 * count the number of occurances for each unique value
 * and sort from most number to least
*/
template <typename T>
inline std::vector<std::pair<T, int>> 
sortNumberOccurances(const std::vector<T> &values) {
	std::unordered_map<T, int> counts;
	for (const auto &i : values) {
		counts[i]++;
	}

	return sortNumberOccurances(counts);
}

template <typename T>
inline std::vector<std::pair<T, int>>
sortNumberOccurances(const std::unordered_map<T, int> &counts) {
	std::vector<std::pair<T, int>> out;
	for (const auto &i : counts) {
		out.push_back(i);
	}

	std::sort(out.begin(), out.end(),
		[](const auto &p1, const auto p2) {
			return (p1.second > p2.second);
		}
	);

	return out;	
}

template <typename T>
inline std::unordered_map<T, int>
countUnique(const std::vector<T> &vec) {
	std::unordered_map<T, int> counts;
	for (const auto &it : vec) {
		counts[it]++;
	}
	return counts;
}

inline std::pair<std::string, std::string>
parseSpace(const std::string &s) {
	// parse a string until delim, and return the parsed string and the rest of the string
	for (int i = 0; i < s.size(); i++) {
		if (isspace(s[i])) {
			return {s.substr(0, i), s.substr(i+1, std::string::npos)};
		}
	}
	return {};
}

inline std::pair<std::string, std::string>
parseDelim(const std::string &s, const char &delim) {
	for (size_t i = 0; i < s.size(); i++) {
		if (s[i] == delim) {
			return {s.substr(0, i), s.substr(i+1, std::string::npos)};
		}
	}
	return {};
}

inline std::string get_extdata(std::string file) {
	Rcpp::Function sys_file("system.file");
	std::string res = Rcpp::as<std::string> (sys_file("extdata", file, Rcpp::_["package"] = "FLAMES"));
	return std::string(res);
}
inline std::string get_tempfile(std::string ext) {
	Rcpp::Function temp_file("tempfile");
	std::string res = Rcpp::as<std::string> (temp_file(Rcpp::_["fileext"] = ext));
	return std::string(res);
}
inline std::string get_tempdir() {
	Rcpp::Function temp_dir("tempdir");
	return Rcpp::as<std::string> (temp_dir());
}
inline void rename_file(std::string from, std::string to) {
	Rcpp::Function rename("file.rename");
	rename(Rcpp::_["from"]=from, Rcpp::_["to"]=to);
}
// inline bool file_exists(std::string file) {
// 	Rcpp::Function fileExists("file.exists");
// 	return Rcpp::as<bool> (fileExists(file));
// }

/*
    checks whether an int vector is strictly increasing
*/
inline bool isStrictlyIncreasing(const std::vector<int> &vec)
{
    bool strictlyIncreasing = 1;
    for (int i = 1; i < (int)vec.size(); i++) {
        if (vec[i] < vec[i - 1]) {
            // then the vector is not strictly increasing
            strictlyIncreasing = 0;
            break;
        }
    }
    return strictlyIncreasing;
}

#endif // UTILITY_H
