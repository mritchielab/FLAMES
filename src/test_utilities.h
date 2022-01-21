#ifndef TEST_UTILITIES_H
#define TEST_UTILITIES_H

#include <string>
#include <vector>
#include <algorithm>

#include <Rcpp.h>

// R system.file
// https://stackoverflow.com/questions/52979101/acessing-data-on-inst-extdata-from-rcpp-catch-tests

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

// compare two unordered streams of values
template <typename T>
inline bool compare_unordered(std::vector<T> a, std::vector<T> b) {
	if (a.size() != b.size()) return false;

	for (int i = 0; i < a.size(); i++) {
		if (std::find(b.begin(), b.end(), a[i]) == b.end()) return false;
	}
	
	return true;
}

// compare two unordered maps
template <typename T, typename U>
inline bool compare_map(std::unordered_map<T, std::vector<U>> a, std::unordered_map<T, std::vector<U>> b) {
	for (const auto &[key, value] : a) {
		if (!compare_unordered<U>(value, b[key])) return false;
	}
	return true;
}
template <typename T, typename U>
inline bool compare_map(std::map<T, std::vector<U>> a, std::map<T, std::vector<U>> b) {
	for (const auto &[key, value] : a) {
		if (!compare_unordered<U>(value, b[key])) return false;
	}
	return true;
}

// Compare two iterable streams, expecting values to be in the same order
template <class Iter>
inline bool compare_stream(Iter &a, Iter &b) {
	for (auto i = a.begin(), j = b.begin(); i != a.end(); i++, j++) {
		if (*i != *j) return false;
	}
	return true;
}




#endif // TEST_UTILITIES_H