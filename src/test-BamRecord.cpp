#include <testthat.h>

#include "utility/bam.h"
#include "test_utilities.h"

#include "classes/BamRecord.h"

// static int testFetch1(const bam1_t *b, void *data) {
// 	int *idx = (int *)data; // grab to int pointer from the void pointer

// 	Rcpp::Rcout << "Qname: " << bam1_qname(b) << "\n";

// 	return 1;
// }

context("BamRecord & utility BAM function testing") {
	// test_that("Generated cigar pairs are the same as pysam") {
	// 	// bam = pysam.AlignmentFile("/Users/voogd.o/Downloads/Test1-ready.bam");
	// 	// reg = [x for x in bam.fetch("chrM", 7000, 7008)][0]
		
	// 	const char *bam_in = "/Users/voogd.o/Downloads/Test1-ready.bam";
	// 	bamFile bam = bam_open(bam_in, "r");
	// 	bam_index_t *bam_index = bam_index_load(bam_in);
	// 	int tid = 0;

	// 	int idx = 0;
	// 	bam_fetch(bam, bam_index, tid, 7000, 7008, &idx, testFetch1);
	// 	bam_close(bam);

	// }
}

// [[Rcpp::export]]
void what2() {
		// const char *bam_in = get_extdata("Test1-ready.bam").c_str();
		// bamFile bam = bam_open(bam_in, "r");
		// bam_index_t *bam_index = bam_index_load(bam_in);
		// int tid = 0;

		// int idx = 0;
		// bam_fetch(bam, bam_index, tid, 7000, 7008, &idx, testFetch1);
		// bam_close(bam);

}