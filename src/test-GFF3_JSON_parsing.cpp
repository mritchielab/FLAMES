#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <functional>
#include <Rcpp.h>
#include <testthat.h>

// #include "ParseGFF3.hpp"
#include "classes/Parser.h"
#include "classes/Pos.h"
#include "classes/StartEndPair.h"
// #include "parse_gene_anno_native.h"
#include "classes/GeneAnnoParser/GeneAnnoParser.h"
#include "classes/Config.h"
#include "classes/GeneAnnoParser/GFFRecord.h"
#include "file-handling/parse_json_config.h"
#include "test_utilities.h"

template <typename T>
bool unorderedMatch(std::vector<T> v1, std::vector<T> v2) {
    for (auto i : v1) {
        bool match = false;
        for (auto j : v2) {
            if (i == j) {
                match = true;
                break;
            }
        }
        if (!match) return false;
    }
    return true;
}

bool compareRecord(GFFRecord &r1, GFFRecord r2) {
	return 
		r1.seqname == r2.seqname && 
		r1.source == r2.source && 
		r1.feature == r2.feature &&
		r1.start == r2.start &&
		r1.end == r2.end &&
		r1.score == r2.score &&
		r1.strand == r2.strand &&
		r1.frame == r2.frame;
}

context("GFF3 File Parsing") {
	test_that("utility parsers sucessfully parse small units") {
		std::string full = "SIRV5	.	gene	1009	11866	.	+	.	ID=gene:SIRV5;gene_id=SIRV5;support_count=129";

		auto p1 = parseSpaces("     \tword");
		expect_true(p1.second =="word");

		auto p2 = parseColumn(full);
		expect_true(p2.first == "SIRV5");

		auto p21 = parseColumn("SIRV6\tgene", '\t');
		expect_true(p21.first == "SIRV6");

		auto p3 = parseKeyValue("one=two");
		expect_true(p3.first == "one");
		expect_true(p3.second == "two");

		auto p4 = parseAttribute("\tID=gene:SIRV5;gene_id=SIRV5;support_count=129");
		expect_true(p4.first == "ID=gene:SIRV5");
		expect_true(p4.second == "gene_id=SIRV5;support_count=129");
	
		auto p5 = parseLine("one\ttwo\tthree\tfour");
		expect_true(p5[0] == "one");
		expect_true(p5[3] == "four");

		auto p51 = parseLine("one\ntwo\tthree", '\n');
		expect_true(p51[0] == "one");
		expect_true(p51[1] == "two\tthree");

		auto p6 = parseUntilChar("awordhere\"rest", '"');
		expect_true(p6.first == "awordhere");
		expect_true(p6.second == "rest");
	}

	test_that("parsed GFF3 files give same result as python script") {
		// file = "/Users/voogd.o/Documents/FLAMESData/data/isoform_annotated.gff3"
		// line1 = SIRV5	.	gene	1009	11866	.	+	.	ID=gene:SIRV5;gene_id=SIRV5;support_count=129
		// x = parseGFF3(file)
		// x.next()
		// >>> GFFRecord(seqid='SIRV5', source=None, type='gene', 1009, 11866, score=None, '+', phase=None, attributes={'support_count': '129', 'ID': 'gene:SIRV5', 'gene_id': 'SIRV5'})
		// R >>> GFFRecord(seqid=SIRV5, source=, type=gene, 1009, 11866, score=-1, +, phase=, attributes={ID: gene:SIRV5,gene_id: SIRV5,support_count: 129,}

		GFFParser p (get_extdata("isoform_annotated.gff3"), true);
		// we know that the file is not empty
		expect_true(!p.isEmpty());

		GFFRecord r1 = p.parseNextRecord();
		expect_true(compareRecord(r1, GFFRecord("SIRV5\t.\tgene\t1001\t11866\t-1\t+\t.\t.", true))); 
		expect_true(r1.attributes["ID"] == "gene:SIRV5"); 
		expect_true(r1.attributes["support_count"] == "173"); 
		expect_true(r1.attributes["gene_id"] == "SIRV5");

		GFFRecord r2 = p.parseNextRecord();
		expect_true(compareRecord(r2, GFFRecord("SIRV5\tnew\ttranscript\t8316\t10991\t-1\t+\t.\t.", true)));
		expect_true(r2.attributes["support_count"] == "12");
		expect_true(r2.attributes["source"] == "new"); 
		expect_true(r2.attributes["transcript_id"] == "SIRV5_8316_10991_1"); 
		expect_true(r2.attributes["ID"] == "transcript:SIRV5_8316_10991_1");
		expect_true(r2.attributes["Parent"] == "gene:SIRV5");
		p.close();
	}
	
	test_that("full gff file parsing yeilds correct results") {
		GeneAnnoParser geneAnnoParser (get_extdata("isoform_annotated.gff3"), true);
		GFFData x = geneAnnoParser.parse();

		// test if sizes of maps are correct
		expect_true(x.chr_to_gene.size() == 7);
		expect_true(x.gene_to_transcript.size() == 7);
		expect_true(x.transcript_dict.size() == 47);
		expect_true(x.transcript_to_exon.size() == 47);


		// test if elements of chr_to_gene are correct
		expect_true(unorderedMatch<std::string>(x.chr_to_gene["SIRV5"], std::vector<std::string> {"SIRV5"}));
		expect_true(unorderedMatch<std::string>(x.chr_to_gene["SIRV4"], std::vector<std::string> {"SIRV4"}));
		expect_true(unorderedMatch<std::string>(x.chr_to_gene["SIRV7"], std::vector<std::string> {"SIRV7"}));
		expect_true(unorderedMatch<std::string>(x.chr_to_gene["SIRV6"], std::vector<std::string> {"SIRV6"}));
		expect_true(unorderedMatch<std::string>(x.chr_to_gene["SIRV1"], std::vector<std::string> {"SIRV1"}));
		expect_true(unorderedMatch<std::string>(x.chr_to_gene["SIRV3"], std::vector<std::string> {"SIRV3"}));
		expect_true(unorderedMatch<std::string>(x.chr_to_gene["SIRV2"], std::vector<std::string> {"SIRV2"}));
		// 
		// test if some elements of transcript_dict are correct
		expect_true(comparePos(x.transcript_dict["SIRV410"], Pos {"SIRV4", 1455, 2771, '+', "SIRV4"}));
		expect_true(comparePos(x.transcript_dict["SIRV606"], Pos {"SIRV6", 2285, 10788, '+', "SIRV6"}));
		expect_true(comparePos(x.transcript_dict["SIRV203"], Pos {"SIRV2", 3665, 5895, '-', "SIRV2"}));
		// test if some elements of gene_to_transcript are correct
		expect_true(
			unorderedMatch<std::string>(
				x.gene_to_transcript["SIRV5"],
				std::vector<std::string> {"SIRV508", "SIRV509", "SIRV503", "SIRV5_8316_10991_1", "SIRV506", "SIRV505", "SIRV512"}
			)
		);
		expect_true(
			unorderedMatch<std::string>(
				x.gene_to_transcript["SIRV6"],
				std::vector<std::string> {"SIRV618", "SIRV606", "SIRV616", "SIRV602", "SIRV604", "SIRV607", "SIRV617", "SIRV609", "SIRV608", "SIRV614", "SIRV615", "SIRV613", "SIRV610", "SIRV611"}
			)
		);

		// test if some elements of transcript_to_exon are correct
		expect_true(
			unorderedMatch<StartEndPair>(
				x.transcript_to_exon["SIRV410"],
				std::vector<StartEndPair> { {1455, 1885}, {2251, 2771} }
			)
		);
		expect_true(
			unorderedMatch<StartEndPair>(
				x.transcript_to_exon["SIRV606"],
				std::vector<StartEndPair> { {2285, 2620}, {2740, 2828}, {3106, 3164}, {10724, 10788} }
			)
		);
		expect_true(
			unorderedMatch<StartEndPair>(
				x.transcript_to_exon["SIRV602"],
				std::vector<StartEndPair> { {1124, 1186}, {1468, 1534}, {1640, 1735}, {2780, 2828}, {3106, 3164}, {10724, 10818}, {11031, 11108}, {11205, 11279} }
			)
		);
	}

	// test_that("full gtf file parses correctly") {
	// 	GeneAnnoParser parser ("/Volumes/voogd.o/FLAMES/inst/extdata/SIRV_anno.gtf", false);
	// 	GFFData x = parser.parse();

	// 	expect_true(x.chr_to_gene.size() == 7);
	// 	expect_true(x.gene_to_transcript.size() == 7);
	// 	expect_true(x.transcript_dict.size() == 69);
	// 	expect_true(x.transcript_to_exon.size() == 69);

	// 	// nned to finish the test cases for this
	// }
}

context("JSON file parsing") {
	test_that("JSON Config parsing builds correct Rcpp::List object and Config object") {
		Rcpp::List c = parse_json_config(get_extdata("SIRV_config_default.json"));
		Config config(c);
		
		PipelineParameters p = config.pipeline_parameters;
		expect_true(p.do_genome_alignment);
		expect_true(p.do_isoform_identification);
		expect_true(p.do_read_realignment);
		expect_true(p.do_transcript_quantification);

		GlobalParameters g = config.global_parameters;
		expect_true(g.generate_raw_isoform);
		expect_true(!g.has_UMI);
		
		IsoformParameters i = config.isoform_parameters;
		expect_true(i.MAX_DIST == 10);
        expect_true(i.MAX_TS_DIST == 100);
        expect_true(i.MAX_SPLICE_MATCH_DIST == 10);
        expect_true(i.MIN_FL_EXON_LEN == 40);
        expect_true(i.MAX_SITE_PER_SPLICE == 3);
        expect_true(i.MIN_SUP_CNT == 10);
        expect_true(i.MIN_CNT_PCT == (float)0.01);
        expect_true(i.MIN_SUP_PCT == (float)0.2);
        expect_true(i.STRAND_SPECIFIC == 1);
        expect_true(i.REMOVE_INCOMP_READS == 5);

		AlignmentParameters a = config.alignment_parameters;
		expect_true(a.use_junctions);
		expect_true(a.no_flank);

		RealignParameters r = config.realign_parameters;
		expect_true(r.use_annotation);

		TranscriptCounting t = config.transcript_counting;
		expect_true(t.min_tr_coverage == 0.75);
		expect_true(t.min_read_coverage == 0.75);
	}
}
