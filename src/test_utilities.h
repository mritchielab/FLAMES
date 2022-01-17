#ifndef TEST_UTILITIES_H
#define TEST_UTILITIES_H

#include <string>
#include <map>
#include <unordered_map>
#include <vector>
#include <set>

#include <Rcpp.h>

// R system.file
// https://stackoverflow.com/questions/52979101/acessing-data-on-inst-extdata-from-rcpp-catch-tests
inline std::string get_extdata(const std::string file) {
	Rcpp::Function sys_file("system.file");
	std::string res = Rcpp::as<std::string> (sys_file("extdata", file, Rcpp::_["package"] = "FLAMES"));
	return std::string(res);
}

inline std::string get_tempfile(const std::string ext) {
	Rcpp::Function temp_file("tempfile");
	std::string res = Rcpp::as<std::string> (temp_file(Rcpp::_["fileext"] = ext));
	return std::string(res);
}

template <typename T>
inline bool compare_vector(const std::vector<T> a, const std::vector<T> b) {
	if (a.size() != b.size()) return false;

	for (int i = 0; i < a.size(); i++) {
		if (a.at(i) != b.at(i)) {
			return false;
		}
	}
	return true;
}

// requires values in each stream to be in the same order
template <class Iter>
inline bool compare_stream(const Iter a, const Iter b) {
	if (a.size() != b.size()) return false;

	for (auto i = a.begin(), j = b.begin(); i != a.end(); i++, j++) {
		if (*i != *j) {
			return false;
		}
	}

	return true;
}

// doesn't expect set values to be in the same orders
// template <typename T>
// inline bool compare_set(const std::set<T> a, const std::set<T> b) {
// 	std::unordered_map<T, int> b_map;
// 	for (auto j = b.begin(); j != b.end(); j++) {
// 		b_map[*j] = 1;
// 	}

// 	for (auto i = a.begin(); i != a.end(); i++) {
// 		if (b_map[*i] != 1) return false;
// 	}

// 	return true;
// }

template <typename T, typename U>
inline bool compare_map(const std::map<T, std::vector<U>> a, const std::map<T, std::vector<U>> b) {
	for (const auto & [key, value] : a) {
		if (!b.count(key)) return false;

		if (!compare_vector(value, b.at(key))) return false;
	}
	return true;
}
template <typename T, typename U>
inline bool compare_map(const std::unordered_map<T, std::vector<U>> a, const std::unordered_map<T, std::vector<U>> b) {
	for (const auto & [key, value] : a) {
		if (!b.count(key)) return false;

		if (!compare_vector(value, b.at(key))) return false;
	}
	return true;
}
#endif // TEST_UTILITIES_H