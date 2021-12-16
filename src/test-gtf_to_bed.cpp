#include <string>
#include <fstream>

#include <testthat.h>

#include "test_utilities.h"
#include "gtf_to_bed.h"
#include "ParseGFF3.hpp"
#include "Parser.h"


context("GTF To Bed file conversion") {
	test_that("full gtf_to_bed works on small file") {
		std::string input = get_extdata("SIRV_anno_test.gtf");
		std::string out = get_tempfile(".bed");

		gtf_to_bed_cpp(input, out, std::string());

		std::ifstream res;
		std::string line;
		res.open(out);
		getline(res, line);
		res.close();
		std::vector<std::string> real = {
			"SIRV1",
			"1000",
			"7814",		
			"SIRV101",
			"1000",
			"-",
			"7552",
			"7814",
			"255,0,0",
			"4",			
			"484,136,253,262",
			"0,5337,5560,6552"
		};
		std::vector<std::string> test = parseLine(line);
		for (int i = 0; i < real.size(); i++) {
			expect_true(real[i] == test[i]);
		}
	}

	test_that("full gtf_to_bed works on medium sample file") {

	}

}

// [[Rcpp::export]]
void tester() {
	std::string input = get_extdata("SIRV_anno_test.gtf");
	std::string out = get_tempfile(".bed");

	gtf_to_bed_cpp(input, out, std::string());

	std::ifstream res;
	std::string line;
	res.open(out);
	getline(res, line);
	res.close();
	std::vector<std::string> real = {
		"SIRV1",
		"1000",
		"7814",		
		"SIRV101",
		"1000",
		"-",
		"7552",
		"7814",
		"255,0,0",
		"4",			
		"484,136,253,262",
		"0,5337,5560,6552"
	};
	std::vector<std::string> test = parseLine(line);
	for (int i = 0; i < real.size(); i++) {
		if (real[i] != test[i]) {
			Rcpp::Rcout << real[i] << " ==? " << test[i] << "\n";
		}
	}
}
