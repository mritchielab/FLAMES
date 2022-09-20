#include <testthat.h>
#include <vector>
#include <unordered_map>
#include <string>
#include <fstream>

#include "classes/StartEndPair.h"
#include "utility/utility.h"
#include "test_utilities.h"

context("Test misc functions & utilities") {
	test_that("Most common elements are found") {
		std::vector<int> x {10, 15, 20, 30, 45, 50, 50, 15, 10, 15};

		expect_true(mostCommon<int>(x) == 15);

		std::vector<std::string> y {"one", "two", "one", "five", "0", "one"};
		expect_true(mostCommon<std::string>(y) == "one");

		std::vector<std::vector<int>> z {
			{1, 2, 3, 4, 5},
			{1, 5, 6, 7, 8},
			{2, 2, 5, 6, 7},
			{3, 3, 5, 4, 7}
		};

		std::vector<int> z1 = mostCommonEachCell(z, 5);
		std::vector<int> z_real{1,2,5,4,7};
		expect_true(compare_stream(z_real, z1));
	}

	test_that("Most Common StartEndPair values are found") {
		std::vector<StartEndPair> x {
			{1, 5}, {2, 10}, {3, 15}, {1, 6}, {1, 10}, {15, 10}
		};

		StartEndPair resX = mostCommonSEP(x);

		expect_true((resX == StartEndPair{1, 10}));
	}

	test_that("can download a small csv file successfully") {
		std::string dest = get_tempfile(".csv");
		std::string in = "https://raw.githubusercontent.com/OliverVoogd/FLAMESData/master/data/small_test.csv";

		// download the file
		download_file(in, dest);

		bool fileExists = file_exists(dest);
		expect_true(fileExists);

		if (fileExists) {
			std::string res[6] {
				"AAGTC",
				"AGTCC",
				"AGTCA",
				"GTACA",
				"AGTAC",
				"AGTCA"
			};

			std::string lines[6];
			std::ifstream file(dest);

			std::string line;
			int i = 0;
			while (getline(file, line)) {
				lines[i++] = line;
			}

			expect_true(i == 6);

			for (int j = 0; j < 6; j++) {
				expect_true(res[j] == lines[j]);
			}
		}
	}
	// TODO: test the rest of the functions in utility/utility.h
}



#include "file-handling/write_tr_to_csv.h"
#include "classes/Pos.h"
#include "file-handling/parse_realigned_bam.h"
// general purpose testing func
// [[Rcpp::export]]
void test1() {
	// std::unordered_map<std::string, int> real_tr_count = {
	// 	{"SIRV608", 30},
	// 	{"SIRV701", 17},
	// 	{"SIRV1_6559_10366_1", 9},
	// 	{"SIRV617", 20},
	// 	{"SIRV615", 43},
	// 	{"SIRV308", 69},
	// 	{"SIRV601", 31},
	// 	{"SIRV201", 16},
	// 	{"SIRV307", 42},
	// 	{"SIRV704", 36},
	// 	{"SIRV611", 21},
	// 	{"SIRV610", 32},
	// 	{"SIRV409", 20},
	// 	{"SIRV206", 18},
	// 	{"SIRV612", 68},
	// 	{"SIRV204", 43},
	// 	{"SIRV502", 13},
	// 	{"SIRV604", 58},
	// 	{"SIRV609", 84},
	// 	{"SIRV614", 29},
	// 	{"SIRV605", 39},
	// 	{"SIRV1_10648_11606_1", 13},
	// 	{"SIRV203", 33},
	// 	{"SIRV706", 24},
	// 	{"SIRV301", 101},
	// 	{"SIRV4_8324_13828_1", 5},
	// 	{"SIRV5_8202_10991_1", 13},
	// 	{"SIRV501", 13},
	// 	{"SIRV607", 42},
	// 	{"SIRV109", 25},
	// 	{"SIRV405", 19},
	// 	{"SIRV410", 142},
	// 	{"SIRV107", 93},
	// 	{"SIRV505", 64},
	// 	{"SIRV404", 48},
	// 	{"SIRV602", 145},
	// 	{"SIRV202", 32},
	// 	{"SIRV507", 20},
	// 	{"SIRV708", 35},
	// 	{"SIRV103", 25},
	// 	{"SIRV616", 140},
	// 	{"SIRV205", 8},
	// 	{"SIRV503", 28},
	// 	{"SIRV606", 37},
	// 	{"SIRV508", 62},
	// 	{"SIRV506", 63},
	// 	{"SIRV102", 16},
	// 	{"SIRV106", 44},
	// 	{"SIRV303", 22},
	// 	{"SIRV509", 27},
	// 	{"SIRV613", 83},
	// 	{"SIRV305", 26},
	// 	{"SIRV302", 52},
	// 	{"SIRV403", 16},
	// 	{"SIRV105", 117},
	// 	{"SIRV408", 20},
	// 	{"SIRV406", 165},
	// 	{"SIRV618", 5},
	// 	{"SIRV108", 27},
	// 	{"SIRV7_56032_147957_1", 11},
	// 	{"SIRV703", 8}
	// };


	std::string realign_bam = "/Users/voogd.o/Documents/FLAMESintermediate/FLAMESP/testfiles/realign2transcript.bam";
	std::string transcript_fa_idx = "/Users/voogd.o/Documents/FLAMESintermediate/FLAMESP/testfiles/transcript_assembly.fa.fai";
	std::string bc_file = "/Users/voogd.o/Documents/FLAMESintermediate/FLAMESP/testfiles/pseudo_barcode_annotation.csv";

	std::unordered_map<std::string, std::unordered_map<std::string, std::vector<std::string>>> bc_tr_count_dict_real = {
		{"sample1.fastq.gz", {
			{"SIRV103", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV106", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV105", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV1_6559_10366_1", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV102", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV107", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV108", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV1_10648_11606_1", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV202", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV201", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV203", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV204", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV206", {"NNN", "NNN", "NNN", "NNN"}},
			{"SIRV301", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV303", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV305", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV3_4004_8940_1", {"NNN", "NNN"}},
			{"SIRV307", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV302", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV3_7659_8940_1", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV308", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV405", {"NNN"}},
			{"SIRV404", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV4_8323_13828_1", {"NNN", "NNN", "NNN"}},
			{"SIRV406", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV409", {"NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV410", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV403", {"NNN", "NNN"}},
			{"SIRV408", {"NNN", "NNN"}},
			{"SIRV505", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV503", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV5_8202_10991_1", {"NNN", "NNN"}},
			{"SIRV508", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV501", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV509", {"NNN", "NNN", "NNN"}},
			{"SIRV511", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV502", {"NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV602", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV611", {"NNN", "NNN"}},
			{"SIRV608", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV604", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV612", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV607", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV609", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV605", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV606", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV616", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV601", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV614", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV610", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV613", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV615", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV617", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV701", {"NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV706", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV704", {"NNN", "NNN", "NNN"}},
			{"SIRV708", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV7_56032_147957_1", {"NNN", "NNN", "NNN", "NNN"}}
		}},
		{"sample2.fastq.gz", {
			{"SIRV103", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV106", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV105", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV102", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV109", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV107", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV108", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV1_10648_11606_1", {"NNN", "NNN", "NNN", "NNN"}},
			{"SIRV202", {"NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV201", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV203", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV205", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV204", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV206", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV301", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV303", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV305", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV307", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV302", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV3_7659_8940_1", {"NNN", "NNN", "NNN", "NNN"}},
			{"SIRV308", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV3_4004_8940_1", {"NNN", "NNN"}},
			{"SIRV405", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV403", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV408", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV404", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV4_8323_13828_1", {"NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV406", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV409", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV410", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV505", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV503", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV5_8202_10991_1", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV508", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV507", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV501", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV509", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV511", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV502", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV602", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV611", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV609", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV608", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV614", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV601", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV612", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV604", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV605", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV607", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV616", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV606", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV610", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV613", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV615", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV618", {"NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV617", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV701", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV703", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV704", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV7_56032_147957_1", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV708", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
			{"SIRV706", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}}
		}}
	};

	// {bc_tr_count_dict, bc_tr_badcov_count_dict, tr_kept}
	RealignedBamData parse_realign = parse_realigned_bam(
		realign_bam,
   		transcript_fa_idx,
    	10,
    	0.75,
   		0.75,
    	{}
	);
	// expected printing of above: 'counted_reads': 2571, 'not_enough_coverage': 417, 'no_good_match': 106

	std::unordered_map<std::string, Pos> transcript_dict_i = {
		{"SIRV1_6559_10366_1", {"SIRV1", 6558, 10366, '-', "SIRV1"}},
		{"SIRV1_10648_11606_1", {"SIRV1", 10647, 11606, '-', "SIRV1"}},
		{"SIRV106", {"SIRV1", 1000, 10786, '-', "SIRV1"}},
		{"SIRV103", {"SIRV1", 1000, 10791, '-', "SIRV1"}},
		{"SIRV108", {"SIRV1", 10582, 11606, '+', "SIRV1"}},
		{"SIRV107", {"SIRV1", 10647, 11643, '-', "SIRV1"}},
		{"SIRV109", {"SIRV1", 10711, 11643, '+', "SIRV1"}},
		{"SIRV105", {"SIRV1", 6449, 10640, '-', "SIRV1"}},
		{"SIRV205", {"SIRV2", 1108, 1631, '+', "SIRV2"}},
		{"SIRV206", {"SIRV2", 4033, 4457, '+', "SIRV2"}},
		{"SIRV202", {"SIRV2", 1035, 5911, '-', "SIRV2"}},
		{"SIRV204", {"SIRV2", 3643, 4732, '-', "SIRV2"}},
		{"SIRV203", {"SIRV2", 3665, 5895, '-', "SIRV2"}},
		{"SIRV3_4004_8940_1", {"SIRV3", 4003, 8940, '+', "SIRV3"}},
		{"SIRV3_7659_8940_1", {"SIRV3", 7658, 8940, '+', "SIRV3"}},
		{"SIRV308", {"SIRV3", 1000, 1982, '-', "SIRV3"}},
		{"SIRV305", {"SIRV3", 4003, 6718, '+', "SIRV3"}},
		{"SIRV301", {"SIRV3", 1944, 8939, '+', "SIRV3"}},
		{"SIRV303", {"SIRV3", 1963, 7822, '+', "SIRV3"}},
		{"SIRV4_8323_13822_1", {"SIRV4", 8322, 13822, '-', "SIRV4"}},
		{"SIRV409", {"SIRV4", 1000, 3403, '+', "SIRV4"}},
		{"SIRV410", {"SIRV4", 1455, 2771, '+', "SIRV4"}},
		{"SIRV406", {"SIRV4", 3637, 5158, '-', "SIRV4"}},
		{"SIRV405", {"SIRV4", 8629, 13937, '-', "SIRV4"}},
		{"SIRV5_8202_10991_1", {"SIRV5", 8201, 10991, '+', "SIRV5"}},
		{"SIRV512", {"SIRV5", 2177, 2406, '-', "SIRV5"}},
		{"SIRV506", {"SIRV5", 1008, 2398, '+', "SIRV5"}},
		{"SIRV508", {"SIRV5", 1008, 10991, '+', "SIRV5"}},
		{"SIRV509", {"SIRV5", 8315, 11866, '+', "SIRV5"}},
		{"SIRV503", {"SIRV5", 8201, 11142, '+', "SIRV5"}},
		{"SIRV617", {"SIRV6", 1544, 1820, '-', "SIRV6"}},
		{"SIRV618", {"SIRV6", 2358, 2547, '-', "SIRV6"}},
		{"SIRV604", {"SIRV6", 1087, 11837, '+', "SIRV6"}},
		{"SIRV602", {"SIRV6", 1124, 11279, '+', "SIRV6"}},
		{"SIRV607", {"SIRV6", 1130, 2540, '+', "SIRV6"}},
		{"SIRV609", {"SIRV6", 1137, 2120, '+', "SIRV6"}},
		{"SIRV611", {"SIRV6", 1303, 1950, '+', "SIRV6"}},
		{"SIRV616", {"SIRV6", 2285, 10788, '+', "SIRV6"}},
		{"SIRV606", {"SIRV6", 2285, 10788, '+', "SIRV6"}},
		{"SIRV610", {"SIRV6", 2472, 11690, '+', "SIRV6"}},
		{"SIRV614", {"SIRV6", 2516, 10815, '+', "SIRV6"}},
		{"SIRV608", {"SIRV6", 3023, 11270, '+', "SIRV6"}},
		{"SIRV613", {"SIRV6", 3105, 11824, '+', "SIRV6"}},
		{"SIRV615", {"SIRV6", 10237, 11330, '+', "SIRV6"}},
		{"SIRV7_56032_147957_1", {"SIRV7", 56031, 147957, '-', "SIRV7"}},
		{"SIRV701", {"SIRV7", 1003, 147923, '-', "SIRV7"}},
		{"SIRV704", {"SIRV7", 55849, 114738, '-', "SIRV7"}},
		{"SIRV706", {"SIRV7", 56031, 147957, '-', "SIRV7"}}
	};
	std::string csv_f = "/Users/voogd.o/Documents/FLAMESintermediate/FLAMESC/transcript_count.csv.gz";
	std::unordered_map<std::string, Pos> transcript_dict_ref = {
		{"SIRV101", {"SIRV1", 1000, 10786, '-', "SIRV1"}},
		{"SIRV102", {"SIRV1", 1006, 10366, '-', "SIRV1"}},
		{"SIRV103", {"SIRV1", 1000, 10791, '-', "SIRV1"}},
		{"SIRV105", {"SIRV1", 6449, 10640, '-', "SIRV1"}},
		{"SIRV106", {"SIRV1", 1000, 10786, '-', "SIRV1"}},
		{"SIRV107", {"SIRV1", 10647, 11643, '-', "SIRV1"}},
		{"SIRV108", {"SIRV1", 10582, 11606, '+', "SIRV1"}},
		{"SIRV109", {"SIRV1", 10711, 11643, '+', "SIRV1"}},
		{"SIRV201", {"SIRV2", 1000, 5907, '-', "SIRV2"}},
		{"SIRV202", {"SIRV2", 1035, 5911, '-', "SIRV2"}},
		{"SIRV203", {"SIRV2", 3665, 5895, '-', "SIRV2"}},
		{"SIRV204", {"SIRV2", 3643, 4732, '-', "SIRV2"}},
		{"SIRV205", {"SIRV2", 1108, 1631, '+', "SIRV2"}},
		{"SIRV206", {"SIRV2", 4033, 4457, '+', "SIRV2"}},
		{"SIRV301", {"SIRV3", 1944, 8939, '+', "SIRV3"}},
		{"SIRV302", {"SIRV3", 1963, 7822, '+', "SIRV3"}},
		{"SIRV303", {"SIRV3", 1963, 7822, '+', "SIRV3"}},
		{"SIRV304", {"SIRV3", 1963, 8937, '+', "SIRV3"}},
		{"SIRV305", {"SIRV3", 4003, 6718, '+', "SIRV3"}},
		{"SIRV306", {"SIRV3", 1944, 8292, '+', "SIRV3"}},
		{"SIRV307", {"SIRV3", 1963, 8939, '+', "SIRV3"}},
		{"SIRV308", {"SIRV3", 1000, 1982, '-', "SIRV3"}},
		{"SIRV309", {"SIRV3", 8797, 9943, '-', "SIRV3"}},
		{"SIRV310", {"SIRV3", 8759, 9914, '-', "SIRV3"}},
		{"SIRV311", {"SIRV3", 4601, 4762, '-', "SIRV3"}},
		{"SIRV403", {"SIRV4", 8322, 15122, '-', "SIRV4"}},
		{"SIRV404", {"SIRV4", 8322, 14623, '-', "SIRV4"}},
		{"SIRV405", {"SIRV4", 8629, 13937, '-', "SIRV4"}},
		{"SIRV406", {"SIRV4", 3637, 5158, '-', "SIRV4"}},
		{"SIRV408", {"SIRV4", 8323, 15122, '-', "SIRV4"}},
		{"SIRV409", {"SIRV4", 1000, 3403, '+', "SIRV4"}},
		{"SIRV410", {"SIRV4", 1455, 2771, '+', "SIRV4"}},
		{"SIRV501", {"SIRV5", 1056, 10991, '+', "SIRV5"}},
		{"SIRV502", {"SIRV5", 1019, 10989, '+', "SIRV5"}},
		{"SIRV503", {"SIRV5", 8201, 11142, '+', "SIRV5"}},
		{"SIRV504", {"SIRV5", 11133, 13606, '+', "SIRV5"}},
		{"SIRV505", {"SIRV5", 1000, 10991, '+', "SIRV5"}},
		{"SIRV506", {"SIRV5", 1008, 2398, '+', "SIRV5"}},
		{"SIRV507", {"SIRV5", 1027, 3598, '+', "SIRV5"}},
		{"SIRV508", {"SIRV5", 1008, 10991, '+', "SIRV5"}},
		{"SIRV509", {"SIRV5", 8315, 11866, '+', "SIRV5"}},
		{"SIRV510", {"SIRV5", 1028, 11867, '+', "SIRV5"}},
		{"SIRV511", {"SIRV5", 1008, 2398, '+', "SIRV5"}},
		{"SIRV512", {"SIRV5", 2177, 2406, '-', "SIRV5"}},
		{"SIRV601", {"SIRV6", 1000, 11826, '+', "SIRV6"}},
		{"SIRV602", {"SIRV6", 1124, 11279, '+', "SIRV6"}},
		{"SIRV603", {"SIRV6", 8999, 10968, '+', "SIRV6"}},
		{"SIRV604", {"SIRV6", 1087, 11837, '+', "SIRV6"}},
		{"SIRV605", {"SIRV6", 1130, 11331, '+', "SIRV6"}},
		{"SIRV606", {"SIRV6", 2285, 10788, '+', "SIRV6"}},
		{"SIRV607", {"SIRV6", 1130, 2540, '+', "SIRV6"}},
		{"SIRV608", {"SIRV6", 3023, 11270, '+', "SIRV6"}},
		{"SIRV609", {"SIRV6", 1137, 2120, '+', "SIRV6"}},
		{"SIRV610", {"SIRV6", 2472, 11690, '+', "SIRV6"}},
		{"SIRV611", {"SIRV6", 1303, 1950, '+', "SIRV6"}},
		{"SIRV612", {"SIRV6", 1087, 11825, '+', "SIRV6"}},
		{"SIRV613", {"SIRV6", 3105, 11824, '+', "SIRV6"}},
		{"SIRV614", {"SIRV6", 2516, 10815, '+', "SIRV6"}},
		{"SIRV615", {"SIRV6", 10237, 11330, '+', "SIRV6"}},
		{"SIRV616", {"SIRV6", 2285, 10788, '+', "SIRV6"}},
		{"SIRV617", {"SIRV6", 1544, 1820, '-', "SIRV6"}},
		{"SIRV618", {"SIRV6", 2358, 2547, '-', "SIRV6"}},
		{"SIRV701", {"SIRV7", 1003, 147923, '-', "SIRV7"}},
		{"SIRV702", {"SIRV7", 1000, 114916, '-', "SIRV7"}},
		{"SIRV703", {"SIRV7", 1000, 147918, '-', "SIRV7"}},
		{"SIRV704", {"SIRV7", 55849, 114738, '-', "SIRV7"}},
		{"SIRV705", {"SIRV7", 1005, 147925, '-', "SIRV7"}},
		{"SIRV706", {"SIRV7", 56031, 147957, '-', "SIRV7"}},
		{"SIRV708", {"SIRV7", 56037, 147957, '-', "SIRV7"}}
	};
	bool has_UMI=false;

	std::unordered_map<std::string, int> tr_count = write_tr_to_csv_cpp(
		parse_realign.bc_tr_count_dict, transcript_dict_i, csv_f, transcript_dict_ref, has_UMI
	);

	Rcpp::Rcout << tr_count.size() << " isoforms in count matrix\n";
}
