#include <string>
#include <vector>
#include <Rcpp.h>

#include "test_utilities.h"

// R system.file
// https://stackoverflow.com/questions/52979101/acessing-data-on-inst-extdata-from-rcpp-catch-tests
std::string get_extdata(std::string file) {
	Rcpp::Function sys_file("system.file");
	std::string res = Rcpp::as<std::string> (sys_file("extdata", file, Rcpp::_["package"] = "FLAMES"));
	return std::string(res);
}

std::string get_tempfile(std::string ext) {
	Rcpp::Function temp_file("tempfile");
	std::string res = Rcpp::as<std::string> (temp_file(Rcpp::_["fileext"] = ext));
	return std::string(res);
}