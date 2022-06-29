#include <testthat.h>

#include "main-functions/find_isoform.h"
#include "file-handling/parse_json_config.h"
#include "test_utilities.h"

// context("Find_isoform main FLAMES process") {
// 	test_that("find_isoform correctly terminates and produces correct BAM file from inputs") {

// 	}
// }

// [[Rcpp::export]]
Rcpp::List test_fi() { // rename to test_find_isoform
	// replace FLAMESData files with downloads?
	std::string gff3 = download_data_file("SIRV_isoforms_multi-fasta-annotation_C_170612a.gtf");
	std::string genome_bam = download_align2genome();
	if (genome_bam == "") {
		return Rcpp::List::create();
	}
    std::string isoform_gff3 = get_tempfile(".gff3");
    std::string tss_tes_stat = get_tempfile(".bedgraph");
	std::string genomefa = download_data_file("SIRV_isoforms_multi-fasta_170612a.fasta");
    std::string transcript_fa = get_tempfile(".fa");
    int         downsample_ratio = 1; 
    Rcpp::List  config_list = parse_json_config(get_extdata("SIRV_config_default.json"));
    std::string raw_splice_isoform = get_tempfile(".gff3");

	// args.gff3 = "/Users/voogd.o/Documents/FlamesNew/FLAMES/inst/extdata/SIRV_anno.gtf"
    // args.fq_dir = "/Users/voogd.o/Documents/FlamesNew/FLAMES/inst/extdata/fastq"
    // args.outdir = "/Users/voogd.o/Documents/FLAMESintermediate/beforeFIPython"
    // args.genomefa = "/Users/voogd.o/Documents/FlamesNew/FLAMES/inst/extdata/SIRV_genomefa.fasta"
    // args.minimap2_dir = "/Users/voogd.o/Documents/GitHub/minimap2"
    // args.config_file = "/Users/voogd.o/Documents/FlamesNew/FLAMES/inst/extdata/SIRV_config_default.json"
    // args.downsample_ratio = 1
    // args.inbam = ""


	Rcpp::List iso = find_isoform(gff3, genome_bam, isoform_gff3, tss_tes_stat, genomefa, 
								   transcript_fa, downsample_ratio, config_list, raw_splice_isoform);

	Rcpp::List res = Rcpp::List::create();
	res.push_back(iso, "isoform");
	res.push_back(isoform_gff3, "isoform_gff3");
	res.push_back(tss_tes_stat, "tss_tes");
	res.push_back(transcript_fa, "transcript_fa");
	res.push_back(raw_splice_isoform, "raw_splice");

	return res;
}