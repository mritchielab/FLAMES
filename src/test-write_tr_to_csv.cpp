#include <vector>
#include <string>
#include <utility>
#include <iostream>
#include <fstream>

#include <testthat.h>
#include <Rcpp.h>

#include "classes/Pos.h"
#include "test_utilities.h"
#include "file-handling/write_tr_to_csv.h"

context("Transcript writing and utility functions") {
	test_that("edit_distance calculates correctly") {
		// add more test cases if these cases aren"t all encompassing
		std::vector<std::pair<std::string, std::string>> ps {
			{"AATGC", "AAGGC"},
			{"AAGGGGC", "AATGA"},
			{"GTCGAGATCGA", "AGATCGTAGC"},
			{"GGGGGGGTGGGCG", "GGGGGGGGGGGAAG"},
			{"AAAAA", "AAAAA"}
		};
		std::vector<int> res = {
			1,
			4,
			7,
			3,
			0
		};

	 	// run the test cases
		bool allgood = true;
		for (int i = 0; i < res.size(); i++) {
			allgood &= edit_distance(ps[i].first, ps[i].second) == res[i];
		}
		expect_true(allgood);
	}

	test_that("umi_dedup dedups correctly") {
		// default case of no UMI
		expect_true(umi_dedup(std::vector<std::string> {"AATGC", "AAGTC"}, false) == 2);

		std::vector<std::string> l {
			"AATGC", "AAGTC", "AACGG", "AGTGC", "GATGC", "AGCAT"
		};
		expect_true(umi_dedup(l, true) == 4);
		// please add more test cases here for more UMIs
	}

	test_that("write_tr_to_csv produces correct csv and byproducts") {
		std::unordered_map<std::string, Pos> transcript_dict_i {
			{"SIRV410", {"SIRV4", 1455, 2771, '+', "SIRV4"}},
			{"SIRV606", {"SIRV6", 2285, 10788, '+', "SIRV6"}},
			{"SIRV5_8202_10991_1", {"SIRV5", 8201, 10991, '+', "SIRV5"}},
			{"SIRV5_8316_10991_1", {"SIRV5", 8315, 10991, '+', "SIRV5"}},
			{"SIRV512", {"SIRV5", 2177, 2406, '-', "SIRV5"}},
			{"SIRV602", {"SIRV6", 1124, 11279, '+', "SIRV6"}},
			{"SIRV604", {"SIRV6", 1087, 11837, '+', "SIRV6"}},
			{"SIRV608", {"SIRV6", 3023, 11270, '+', "SIRV6"}},
			{"SIRV609", {"SIRV6", 1137, 2120, '+', "SIRV6"}},
			{"SIRV701", {"SIRV7", 1003, 147923, '-', "SIRV7"}},
			{"SIRV704", {"SIRV7", 55849, 114738, '-', "SIRV7"}},
			{"SIRV706", {"SIRV7", 56031, 147957, '-', "SIRV7"}},
			{"SIRV4_8323_13822_1", {"SIRV4", 8322, 13822, '-', "SIRV4"}},
			{"SIRV301", {"SIRV3", 1944, 8939, '+', "SIRV3"}},
			{"SIRV7_56032_147957_1", {"SIRV7", 56031, 147957, '-', "SIRV7"}},
			{"SIRV303", {"SIRV3", 1963, 7822, '+', "SIRV3"}},
			{"SIRV305", {"SIRV3", 4003, 6718, '+', "SIRV3"}},
			{"SIRV508", {"SIRV5", 1008, 10991, '+', "SIRV5"}},
			{"SIRV509", {"SIRV5", 8315, 11866, '+', "SIRV5"}},
			{"SIRV406", {"SIRV4", 3637, 5158, '-', "SIRV4"}},
			{"SIRV405", {"SIRV4", 8629, 13937, '-', "SIRV4"}},
			{"SIRV503", {"SIRV5", 8201, 11142, '+', "SIRV5"}},
			{"SIRV409", {"SIRV4", 1000, 3403, '+', "SIRV4"}},
			{"SIRV506", {"SIRV5", 1008, 2398, '+', "SIRV5"}},
			{"SIRV205", {"SIRV2", 1108, 1631, '+', "SIRV2"}},
			{"SIRV204", {"SIRV2", 3643, 4732, '-', "SIRV2"}},
			{"SIRV206", {"SIRV2", 4033, 4457, '+', "SIRV2"}},
			{"SIRV308", {"SIRV3", 1000, 1982, '-', "SIRV3"}},
			{"SIRV203", {"SIRV2", 3665, 5895, '-', "SIRV2"}},
			{"SIRV202", {"SIRV2", 1035, 5911, '-', "SIRV2"}},
			{"SIRV616", {"SIRV6", 2285, 10788, '+', "SIRV6"}},
			{"SIRV617", {"SIRV6", 1544, 1820, '-', "SIRV6"}},
			{"SIRV614", {"SIRV6", 2516, 10815, '+', "SIRV6"}},
			{"SIRV615", {"SIRV6", 10237, 11330, '+', "SIRV6"}},
			{"SIRV613", {"SIRV6", 3105, 11824, '+', "SIRV6"}},
			{"SIRV610", {"SIRV6", 2472, 11690, '+', "SIRV6"}},
			{"SIRV611", {"SIRV6", 1303, 1950, '+', "SIRV6"}},
			{"SIRV108", {"SIRV1", 10582, 11606, '+', "SIRV1"}},
			{"SIRV109", {"SIRV1", 10711, 11643, '+', "SIRV1"}},
			{"SIRV106", {"SIRV1", 1000, 10786, '-', "SIRV1"}},
			{"SIRV107", {"SIRV1", 10647, 11643, '-', "SIRV1"}},
			{"SIRV105", {"SIRV1", 6449, 10640, '-', "SIRV1"}},
			{"SIRV103", {"SIRV1", 1000, 10791, '-', "SIRV1"}}
		};
		std::unordered_map<std::string, Pos> transcript_dict {
			{"SIRV410", {"SIRV4", 1455, 2771, '+', "SIRV4"}}, 
			{"SIRV706", {"SIRV7", 56031, 147957, '-', "SIRV7"}}, 
			{"SIRV618", {"SIRV6", 2358, 2547, '-', "SIRV6"}}, 
			{"SIRV511", {"SIRV5", 1008, 2398, '+', "SIRV5"}}, 
			{"SIRV510", {"SIRV5", 1028, 11867, '+', "SIRV5"}}, 
			{"SIRV512", {"SIRV5", 2177, 2406, '-', "SIRV5"}}, 
			{"SIRV601", {"SIRV6", 1000, 11826, '+', "SIRV6"}}, 
			{"SIRV616", {"SIRV6", 2285, 10788, '+', "SIRV6"}}, 
			{"SIRV603", {"SIRV6", 8999, 10968, '+', "SIRV6"}}, 
			{"SIRV509", {"SIRV5", 8315, 11866, '+', "SIRV5"}}, 
			{"SIRV605", {"SIRV6", 1130, 11331, '+', "SIRV6"}},
			{"SIRV604", {"SIRV6", 1087, 11837, '+', "SIRV6"}},
			{"SIRV607", {"SIRV6", 1130, 2540, '+', "SIRV6"}},
			{"SIRV606", {"SIRV6", 2285, 10788, '+', "SIRV6"}},
			{"SIRV609", {"SIRV6", 1137, 2120, '+', "SIRV6"}},
			{"SIRV608", {"SIRV6", 3023, 11270, '+', "SIRV6"}},
			{"SIRV702", {"SIRV7", 1000, 114916, '-', "SIRV7"}},
			{"SIRV703", {"SIRV7", 1000, 147918, '-', "SIRV7"}},
			{"SIRV704", {"SIRV7", 55849, 114738, '-', "SIRV7"}},
			{"SIRV614", {"SIRV6", 2516, 10815, '+', "SIRV6"}},
			{"SIRV311", {"SIRV3", 4601, 4762, '-', "SIRV3"}},
			{"SIRV310", {"SIRV3", 8759, 9914, '-', "SIRV3"}},
			{"SIRV615", {"SIRV6", 10237, 11330, '+', "SIRV6"}},
			{"SIRV705", {"SIRV7", 1005, 147925, '-', "SIRV7"}},
			{"SIRV701", {"SIRV7", 1003, 147923, '-', "SIRV7"}},
			{"SIRV613", {"SIRV6", 3105, 11824, '+', "SIRV6"}},
			{"SIRV610", {"SIRV6", 2472, 11690, '+', "SIRV6"}},
			{"SIRV611", {"SIRV6", 1303, 1950, '+', "SIRV6"}},
			{"SIRV403", {"SIRV4", 8322, 15122, '-', "SIRV4"}},
			{"SIRV508", {"SIRV5", 1008, 10991, '+', "SIRV5"}},
			{"SIRV501", {"SIRV5", 1056, 10991, '+', "SIRV5"}},
			{"SIRV406", {"SIRV4", 3637, 5158, '-', "SIRV4"}},
			{"SIRV405", {"SIRV4", 8629, 13937, '-', "SIRV4"}},
			{"SIRV404", {"SIRV4", 8322, 14623, '-', "SIRV4"}},
			{"SIRV502", {"SIRV5", 1019, 10989, '+', "SIRV5"}},
			{"SIRV503", {"SIRV5", 8201, 11142, '+', "SIRV5"}},
			{"SIRV409", {"SIRV4", 1000, 3403, '+', "SIRV4"}},
			{"SIRV408", {"SIRV4", 8323, 15122, '-', "SIRV4"}},
			{"SIRV506", {"SIRV5", 1008, 2398, '+', "SIRV5"}},
			{"SIRV507", {"SIRV5", 1027, 3598, '+', "SIRV5"}},
			{"SIRV504", {"SIRV5", 11133, 13606, '+', "SIRV5"}},
			{"SIRV505", {"SIRV5", 1000, 10991, '+', "SIRV5"}},
			{"SIRV205", {"SIRV2", 1108, 1631, '+', "SIRV2"}},
			{"SIRV204", {"SIRV2", 3643, 4732, '-', "SIRV2"}},
			{"SIRV617", {"SIRV6", 1544, 1820, '-', "SIRV6"}},
			{"SIRV206", {"SIRV2", 4033, 4457, '+', "SIRV2"}},
			{"SIRV308", {"SIRV3", 1000, 1982, '-', "SIRV3"}},
			{"SIRV309", {"SIRV3", 8797, 9943, '-', "SIRV3"}},
			{"SIRV203", {"SIRV2", 3665, 5895, '-', "SIRV2"}},
			{"SIRV202", {"SIRV2", 1035, 5911, '-', "SIRV2"}},
			{"SIRV304", {"SIRV3", 1963, 8937, '+', "SIRV3"}},
			{"SIRV305", {"SIRV3", 4003, 6718, '+', "SIRV3"}},
			{"SIRV306", {"SIRV3", 1944, 8292, '+', "SIRV3"}},
			{"SIRV307", {"SIRV3", 1963, 8939, '+', "SIRV3"}},
			{"SIRV612", {"SIRV6", 1087, 11825, '+', "SIRV6"}},
			{"SIRV301", {"SIRV3", 1944, 8939, '+', "SIRV3"}},
			{"SIRV302", {"SIRV3", 1963, 7822, '+', "SIRV3"}},
			{"SIRV303", {"SIRV3", 1963, 7822, '+', "SIRV3"}},
			{"SIRV708", {"SIRV7", 56037, 147957, '-', "SIRV7"}},
			{"SIRV108", {"SIRV1", 10582, 11606, '+', "SIRV1"}},
			{"SIRV109", {"SIRV1", 10711, 11643, '+', "SIRV1"}},
			{"SIRV106", {"SIRV1", 1000, 10786, '-', "SIRV1"}},
			{"SIRV107", {"SIRV1", 10647, 11643, '-', "SIRV1"}},
			{"SIRV105", {"SIRV1", 6449, 10640, '-', "SIRV1"}},
			{"SIRV102", {"SIRV1", 1006, 10366, '-', "SIRV1"}},
			{"SIRV103", {"SIRV1", 1000, 10791, '-', "SIRV1"}},
			{"SIRV101", {"SIRV1", 1000, 10786, '-', "SIRV1"}},
			{"SIRV602", {"SIRV6", 1124, 11279, '+', "SIRV6"}},
			{"SIRV201", {"SIRV2", 1000, 5907, '-', "SIRV2"}}
		};
		std::unordered_map<std::string, std::unordered_map<std::string, std::vector<std::string>>> bc_tr_count_dict {
			{"sample2.fastq.gz", {
					{"SIRV410", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}}, 
					{"SIRV618", {"NNN", "NNN", "NNN", "NNN", "NNN"}},
					{"SIRV408", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
					{"SIRV5_8202_10991_1", {"NNN", "NNN", "NNN", "NNN"}},
					{"SIRV5_8316_10991_1", {"NNN", "NNN", "NNN"}},
					{"SIRV708", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
					{"SIRV602", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
					{"SIRV604", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
					{"SIRV607", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
					{"SIRV606", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
					{"SIRV609", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
					{"SIRV608", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
					{"SIRV703", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
					{"SIRV704", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
					{"SIRV706", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
					{"SIRV4_8323_13822_1", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
					{"SIRV615", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
					{"SIRV701", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
					{"SIRV301", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
					{"SIRV7_56032_147957_1", {"NNN"}},
					{"SIRV303", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
					{"SIRV403", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
					{"SIRV508", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
					{"SIRV509", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
					{"SIRV406", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
					{"SIRV405", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
					{"SIRV404", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
					{"SIRV502", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
					{"SIRV503", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
					{"SIRV409", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
					{"SIRV501", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
					{"SIRV506", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
					{"SIRV507", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
					{"SIRV505", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
					{"SIRV205", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
					{"SIRV204", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
					{"SIRV617", {"NNN", "NNN", "NNN", "NNN", "NNN"}},
					{"SIRV206", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
					{"SIRV308", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
					{"SIRV203", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
					{"SIRV202", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
					{"SIRV616", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
					{"SIRV305", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
					{"SIRV614", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
					{"SIRV307", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
					{"SIRV612", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
					{"SIRV613", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
					{"SIRV610", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
					{"SIRV611", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
					{"SIRV108", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
					{"SIRV109", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
					{"SIRV106", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
					{"SIRV107", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
					{"SIRV105", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
					{"SIRV102", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
					{"SIRV103", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
					{"SIRV201", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}}
				}
			},
			{"Sample1.fastq.gz", {
					{"SIRV410", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
					{"SIRV403", {"NNN", "NNN"}},
					{"SIRV408", {"NNN", "NNN"}},
					{"SIRV5_8202_10991_1", {"NNN"}},
					{"SIRV5_8316_10991_1", {"NNN"}},
					{"SIRV708", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
					{"SIRV602", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
					{"SIRV604", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
					{"SIRV607", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
					{"SIRV606", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
					{"SIRV609", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
					{"SIRV701", {"NNN", "NNN", "NNN", "NNN", "NNN"}},
					{"SIRV704", {"NNN", "NNN", "NNN"}},
					{"SIRV706", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
					{"SIRV4_8323_13822_1", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
					{"SIRV615", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
					{"SIRV613", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
					{"SIRV303", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
					{"SIRV608", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
					{"SIRV508", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
					{"SIRV509", {"NNN", "NNN"}},
					{"SIRV406", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
					{"SIRV405", {"NNN"}},
					{"SIRV404", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
					{"SIRV502", {"NNN", "NNN", "NNN", "NNN", "NNN"}},
					{"SIRV503", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
					{"SIRV409", {"NNN", "NNN", "NNN", "NNN", "NNN"}},
					{"SIRV501", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
					{"SIRV506", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
					{"SIRV505", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
					{"SIRV204", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
					{"SIRV617", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
					{"SIRV206", {"NNN", "NNN", "NNN", "NNN"}},
					{"SIRV308", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
					{"SIRV203", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
					{"SIRV202", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
					{"SIRV616", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
					{"SIRV305", {"NNN", "NNN", "NNN", "NNN"}},
					{"SIRV614", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
					{"SIRV307", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
					{"SIRV612", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
					{"SIRV301", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
					{"SIRV610", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
					{"SIRV611", {"NNN", "NNN"}},
					{"SIRV108", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
					{"SIRV106", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
					{"SIRV107", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
					{"SIRV105", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
					{"SIRV102", {"NNN", "NNN", "NNN", "NNN", "NNN"}},
					{"SIRV103", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}},
					{"SIRV201", {"NNN", "NNN", "NNN", "NNN", "NNN", "NNN", "NNN"}}
				}
			}
		};
		
		std::string out_csv = get_tempfile(".csv");
		std::cout << out_csv << "\n";
		std::unordered_map<std::string, int> result {
			{"SIRV410,SIRV4,23,119", 1},
			{"SIRV403,SIRV4,14,2", 1},
			{"SIRV618,SIRV6,5,0", 1},
			{"SIRV408,SIRV4,19,2", 1},
			{"SIRV5_8202_10991_1,SIRV5,4,1", 1},
			{"SIRV5_8316_10991_1,SIRV5,3,1", 1},
			{"SIRV708,SIRV7,15,20", 1},
			{"SIRV602,SIRV6,24,122", 1},
			{"SIRV604,SIRV6,26,34", 1},
			{"SIRV607,SIRV6,14,21", 1},
			{"SIRV606,SIRV6,24,13", 1},
			{"SIRV609,SIRV6,39,38", 1},
			{"SIRV701,SIRV7,12,5", 1},
			{"SIRV703,SIRV7,8,0", 1},
			{"SIRV704,SIRV7,33,3", 1},
			{"SIRV706,SIRV7,13,11", 1},
			{"SIRV4_8323_13822_1,SIRV4,14,6", 1},
			{"SIRV615,SIRV6,36,8", 1},
			{"SIRV613,SIRV6,22,20", 1},
			{"SIRV7_56032_147957_1,SIRV7,1,0", 1},
			{"SIRV611,SIRV6,19,2", 1},
			{"SIRV608,SIRV6,26,7", 1},
			{"SIRV508,SIRV5,12,50", 1},
			{"SIRV509,SIRV5,24,2", 1},
			{"SIRV406,SIRV4,39,126", 1},
			{"SIRV405,SIRV4,18,1", 1},
			{"SIRV404,SIRV4,15,31", 1},
			{"SIRV502,SIRV5,7,5", 1},
			{"SIRV503,SIRV5,19,7", 1},
			{"SIRV409,SIRV4,15,5", 1},
			{"SIRV501,SIRV5,23,9", 1},
			{"SIRV506,SIRV5,38,24", 1},
			{"SIRV507,SIRV5,20,0", 1},
			{"SIRV505,SIRV5,23,41", 1},
			{"SIRV205,SIRV2,8,0", 1},
			{"SIRV204,SIRV2,19,24", 1},
			{"SIRV617,SIRV6,5,6", 1},
			{"SIRV206,SIRV2,14,4", 1},
			{"SIRV308,SIRV3,12,57", 1},
			{"SIRV203,SIRV2,24,6", 1},
			{"SIRV202,SIRV2,16,49", 1},
			{"SIRV616,SIRV6,24,117", 1},
			{"SIRV305,SIRV3,19,4", 1},
			{"SIRV614,SIRV6,20,9", 1},
			{"SIRV307,SIRV3,18,25", 1},
			{"SIRV612,SIRV6,25,64", 1},
			{"SIRV301,SIRV3,14,89", 1},
			{"SIRV610,SIRV6,25,10", 1},
			{"SIRV303,SIRV3,19,14", 1},
			{"SIRV108,SIRV1,20,8", 1},
			{"SIRV109,SIRV1,25,0", 1},
			{"SIRV106,SIRV1,18,24", 1},
			{"SIRV107,SIRV1,9,87", 1},
			{"SIRV105,SIRV1,17,124", 1},
			{"SIRV102,SIRV1,8,5", 1},
			{"SIRV103,SIRV1,7,17", 1},
			{"SIRV201,SIRV2,9,7", 1}
		};

		std::unordered_map<std::string, int> tr_cnt = write_tr_to_csv_cpp(bc_tr_count_dict, transcript_dict_i, out_csv, transcript_dict, false);

		// test the resulting file
		std::ifstream output (out_csv);
		std::string line;
		std::getline(output, line); // ignore the header line
		int numlines = 0;
		bool alllinesgood = true;
		while(std::getline(output, line)) {
			if (line.length()) {
				numlines++;
				alllinesgood = alllinesgood && result.count(line);
			}
		}
		expect_true(alllinesgood);
		expect_true(numlines == result.size());

		// test the retured value
		std::unordered_map<std::string, int> real_cnt {
			{"SIRV410", 142},
			{"SIRV203", 30},
			{"SIRV5_8202_10991_1", 5},
			{"SIRV5_8316_10991_1", 4},
			{"SIRV708", 35},
			{"SIRV602", 146},
			{"SIRV604", 60},
			{"SIRV607", 35},
			{"SIRV606", 37},
			{"SIRV609", 77},
			{"SIRV608", 33},
			{"SIRV703", 8},
			{"SIRV704", 36},
			{"SIRV7_56032_147957_1", 1},
			{"SIRV4_8323_13822_1", 20},
			{"SIRV307", 43},
			{"SIRV701", 17},
			{"SIRV301", 103},
			{"SIRV706", 24},
			{"SIRV303", 33},
			{"SIRV305", 23},
			{"SIRV403", 16},
			{"SIRV508", 62},
			{"SIRV509", 26},
			{"SIRV406", 165},
			{"SIRV405", 19},
			{"SIRV404", 46},
			{"SIRV502", 12},
			{"SIRV503", 26},
			{"SIRV409", 20},
			{"SIRV408", 21},
			{"SIRV506", 62},
			{"SIRV507", 20},
			{"SIRV505", 64},
			{"SIRV205", 8},
			{"SIRV204", 43},
			{"SIRV206", 18},
			{"SIRV308", 69},
			{"SIRV618", 5},
			{"SIRV202", 65},
			{"SIRV616", 141},
			{"SIRV617", 11},
			{"SIRV614", 29},
			{"SIRV615", 44},
			{"SIRV612", 89},
			{"SIRV613", 42},
			{"SIRV610", 35},
			{"SIRV611", 21},
			{"SIRV108", 28},
			{"SIRV109", 25},
			{"SIRV106", 42},
			{"SIRV107", 96},
			{"SIRV105", 141},
			{"SIRV102", 13},
			{"SIRV103", 24},
			{"SIRV501", 32},
			{"SIRV201", 16}
		};
		
		bool allgood = true;
		for (const auto & [key, value] : real_cnt) {
			allgood &= tr_cnt[key] == value;
		}
		expect_true(allgood);
	}
}
