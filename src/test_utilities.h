#ifndef TEST_UTILITIES_H
#define TEST_UTILITIES_H

#include <string>
#include <vector>
#include <algorithm>
#include <functional>

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
inline bool compare_unordered(const std::vector<T> &a, const std::vector<T> &b) {
	if (a.size() != b.size()) return false;

	for (int i = 0; i < a.size(); i++) {
		if (std::find(b.begin(), b.end(), a[i]) == b.end()) return false;
	}
	
	return true;
}

// compare two unordered maps
template <typename T, typename U>
inline bool compare_map(const std::unordered_map<T, std::vector<U>> &a, const std::unordered_map<T, std::vector<U>> &b) {
	for (const auto &[key, value] : a) {
		if (!compare_unordered<U>(value, b.at(key))) return false;
	}
	return true;
}
template <typename T, typename U>
inline bool compare_map(const std::map<T, std::vector<U>> &a, const std::map<T, std::vector<U>> &b) {
	for (const auto &[key, value] : a) {
		if (!compare_unordered<U>(value, b.at(key))) return false;
	}
	return true;
}

// Compare two iterable streams, expecting values to be in the same order
template <class Iter>
inline bool compare_stream(const Iter &a, const Iter &b) {
	for (auto i = a.begin(), j = b.begin(); i != a.end(); i++, j++) {
		if (*i != *j) return false;
	}
	return true;
}

template <typename T, typename U>
inline std::vector<U> map(const std::vector<T> &a, std::function<U(const T &)> f) {
	std::vector<U> out;
	for (auto it : a) {
		out.push_back(f(it));
	}
	return out;
}



#endif // TEST_UTILITIES_H