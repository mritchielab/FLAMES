#ifndef UTILITY_H
#define UTILITY_H

#include <string>
#include <sstream>
#include <unordered_map>
#include <vector>
#include <functional>
#include <algorithm>
#include <set>

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
	
	template <typename T, class Comp>
	void sort(std::vector<T> &vec, Comp f) {
		std::sort(vec.begin(), vec.end(), f);
	}

	// template <typename T, typename U>
	// U reduce(const std::vector<T> &range, std::function<U(const T &, const U &) f) {
	// 	U res = range[0];
	// 	for (int i = 1; i < range.size(); i++) {
	// 		res = f(range[i], res);
	// 	}
	// 	return res;
	// }
	// template <typename T, typename U>
	// U reduce(const std::vector<T> &range, U init, std::function<U(const T &, const U &) f) {
	// 	U res = init;
	// 	for (int i = 0; i < range.size(); i++) {
	// 		res = f(range[i], res);
	// 	}
	// 	return res;
	// }
}


template <typename T>
void printMap(std::unordered_map<std::string, T> &, std::function<std::string (T)>);

std::string PosToStr(Pos);

std::string VecToStr(std::vector<std::string>);

template <typename T>
inline std::string VecToStr(const std::vector<T> &v) {
	std::stringstream ss;
	ss << "{";
	for (const auto &it : v) {
		ss << it << ",";
	}
	ss << "}";
	return ss.str();
}
/*
 * At vector position, determine the most value across all vectors
 * @return a vector of the most common elements at each position
*/
template <typename T>
static inline std::vector<T> 
mostCommonEachCell(const std::vector<std::vector<T>> &values, int maxSize) {
	std::vector<T> newValue;
	for (int i = 0; i < maxSize; i++) {
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
static inline T mostCommon(const std::vector<T> &values) {
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

StartEndPair mostCommonSEP(const std::vector<StartEndPair> &);

/*
 * Taking a vector of values (which contains duplicates)
 * count the number of occurances for each unique value
 * and sort from most number to least
*/
template <typename T>
std::vector<std::pair<T, int>> 
static inline sortNumberOccurances(const std::vector<T> &values) {
	std::unordered_map<T, int> counts;
	for (const auto &i : values) {
		counts[i]++;
	}

	return sortNumberOccurances(counts);
}

template <typename T>
std::vector<std::pair<T, int>>
static inline sortNumberOccurances(const std::unordered_map<T, int> &counts) {
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
std::unordered_map<T, int>
static inline countUnique(const std::vector<T> &vec) {
	std::unordered_map<T, int> counts;
	for (const auto &it : vec) {
		counts[it]++;
	}
	return counts;
}

#endif // UTILITY_H
