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