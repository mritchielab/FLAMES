#include <testthat.h>
#include <vector>
#include <algorithm>

#include "classes/StartEndPair.h"
#include "classes/GeneBlocks.h"

template <typename T>
bool vecEqual(std::vector<T> a, std::vector<T> b) {
	if (a.size() != b.size()) return false;
	for (int i = 0; i < a.size(); i++) {
		if (a[i] != b[i]) {
			return false;
		}
	}
	return true;
}

context("Data Structure tests") {
	// add more tests for any auxillirary DS note tested elsewhere

	test_that("StartEndPair operations allow for correct sorting") {
		std::vector<StartEndPair> v = {
			{3, 4},
			{1, 2},
			{9, 10},
			{5, 6},
			{7, 8}
		};

		expect_true(v[0] < v[2]);
		expect_true(!(v[0] < v[1]));
		expect_true(v[0] <= (StartEndPair {3, 4}));
		expect_true(v[3] > v[1]);
		expect_true(!(v[3] > v[4]));
		expect_true(v[3] >= (StartEndPair {5, 6}));
		expect_true(v[1] == (StartEndPair {1, 2}));

		std::sort(v.begin(), v.end());
		std::vector<StartEndPair> sorted = {
			{1, 2},
			{3, 4},
			{5, 6},
			{7, 8},
			{9, 10}
		};
		expect_true(vecEqual<StartEndPair>(v, sorted));
	}

	test_that("GeneBlocks builds and adds genes correctly") {
		GeneBlocks gb(0, 100, std::vector<std::string> (), "g1");

		expect_true(gb.gene_to_transcript["g1"].size() == 0);

		gb.add_gene(20, 165, std::vector<std::string> {
			"one", "two", "three", "four"
		}, "g2");

		expect_true(gb.end == 165);
		expect_true(gb.start == 0);
		expect_true(gb.transcript_list.size() == 4);
		expect_true(gb.gene_to_transcript["g2"].size() == 4);

		
	}

}
