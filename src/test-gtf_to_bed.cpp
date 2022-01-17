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
		std::string input = get_extdata("SIRV_anno.gtf");
		std::string out = get_tempfile(".bed");

		gtf_to_bed_cpp(input, out, std::string());

		std::ifstream res;
		res.open(out);
		std::string line;
		std::unordered_map<std::string, int> m;
		while (getline(res, line)) {
			m[line] = 100;
		}

		expect_true(m.size() == 69);
		// hand test a few lines
		expect_true(m["SIRV1\t1000\t10786\tSIRV101\t1000\t-\t1000\t10786\t255,0,0\t6\t484,136,253,262,84,342\t0,5337,5560,6552,9282,9444"]);
		expect_true(m["SIRV7\t56037\t147957\tSIRV708\t1000\t-\t147608\t147957\t255,0,0\t6\t60,104,67,35,274,349\t0,14846,22804,22891,58649,91571"]);
		expect_true(m["SIRV6\t2285\t10788\tSIRV616\t1000\t+\t2285\t10788\t255,0,0\t4\t335,74,58,64\t0,455,821,8439"]);
		expect_true(m["SIRV3\t4601\t4762\tSIRV311\t1000\t-\t4601\t4762\t255,0,0\t1\t161\t0"]);
		expect_true(m["SIRV2\t1035\t5911\tSIRV202\t1000\t-\t1035\t5911\t255,0,0\t11\t626,112,91,128,129,220,160,128,141,113,123\t0,706,938,1639,1846,2070,2630,2931,3303,3652,4753"]);
	}

}

