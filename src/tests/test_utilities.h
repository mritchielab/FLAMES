#ifndef TEST_UTILITIES_H
#define TEST_UTILITIES_H

#include <string>
#include <vector>
#include <algorithm>
#include <functional>

#include <Rcpp.h>

#include "../utility/utility.h"
// R system.file
// https://stackoverflow.com/questions/52979101/acessing-data-on-inst-extdata-from-rcpp-catch-tests

inline void download_file(std::string url, std::string dest) {
	Rcpp::Function downloadFile("download.file");
	downloadFile(Rcpp::Named("url")=url, Rcpp::Named("destfile")=dest);
}
inline std::string download_data_file(std::string name) {
	std::string tempfile = get_tempfile(name);
	download_file("https://raw.githubusercontent.com/OliverVoogd/FLAMESData/master/data/" + name, tempfile);
	return tempfile;
}

inline std::string download_align2genome() {
	std::string bam_dir = get_tempdir();
	std::string bam_in = bam_dir + "/align.bam";
	std::string bam_bai = bam_in + ".bai";
	std::string file_url = "https://raw.githubusercontent.com/OliverVoogd/FLAMESData/master/data/align2genome.bam";
	download_file(file_url, bam_in); 
	download_file(file_url + ".bai", bam_bai);
	// std::string bam_in = "/Users/voogd.o/Documents/FLAMESData/data/align2genome.bam";
	if (!file_exists(bam_in) || !file_exists(bam_bai)) {
		Rcpp::Rcout << "required BAM file or BAM index file could not be downloaded from " << file_url << "\n";
		return "";
	}

	return bam_in;
}
// compare two unordered streams of values
template <typename T>
inline bool compare_unordered(const std::vector<T> &a, const std::vector<T> &b) {
	if (a.size() != b.size()) return false;

	for (int i = 0; i < (int)a.size(); i++) {
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
// compare two unordered maps of single values
template <typename T, typename U>
inline bool compare_map(const std::unordered_map<T, U> &a, const std::unordered_map<T, U> &b) {
	for (const auto &[key, value] : a) {
		if (value != b.at(key)) return false;
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

template <typename K, typename T, typename U>
inline std::unordered_map<K, U> map(const std::unordered_map<K, T> &a, std::function<U(const T &)> f) {
	std::unordered_map<K, U> out;
	for (const auto &[key, value] : a) {
		out[key] = f(value);
	}
	return out;
}

template <typename T>
inline std::string printVec(const std::vector<T> &a) {
	std::stringstream ss;
	ss << "{" << *(a.begin());
	for (auto it = a.begin()+1; it != a.end(); it++) {
		ss << ", " << *it;
	}
	ss << "}";
	return ss.str();
}
inline std::string printVecStr(const std::vector<std::string> &a) {
	if (a.size() == 0) return "";
	std::stringstream ss;
	ss << "{" << a[0];
	for (int i = 1; i < (int)a.size(); i++) {
		ss << ", " << a[i];
	}
	ss << "}";
	return ss.str();
}
#endif // TEST_UTILITIES_H