#include <string>
#include <vector>
#include <unordered_map>
#include <fstream>
#include <testthat.h>
#include <Rcpp.h>

#include "test_utilities.h"

#include "main-functions/group_bam2isoform.h"
#include "classes/Config.h"
#include "classes/GFFData.h"
#include "classes/GeneBlocks.h"
#include "classes/GeneAnnoParser/GeneAnnoParser.h"
#include "classes/Pos.h"
#include "classes/StartEndPair.h"
#include "classes/Parser.h"
#include "classes/BamRecord.h"
#include "classes/DataStruct.h"

static int
fetch_function(const bam1_t *b, void *data)
{
    DataStruct * data_struct = (DataStruct*)data;

    std::vector<BAMRecord> * records = data_struct->records;
    bam_header_t * header = data_struct->header;

    BAMRecord rec = read_record(b, header);
    records->push_back(rec);
	return 0;
}

context("Group BAM 2 Isoform") {
	test_that("get_blocks produces correct blocks") {
		// initial test of single position
		// dataset provided by https://github.com/brainstorm/tiny-test-data
		std::string bamIn = get_extdata("Test1-ready.bam");
		bamFile bam = bam_open(bamIn.c_str(), "r");
		bam_index_t *bam_index = bam_index_load(bamIn.c_str());
		bam_header_t *header = bam_header_read(bam);
		int tid = bam_get_tid(header, "chrM");

		// fetch qname ERR032227.10543296
		std::vector<BAMRecord> records;
		DataStruct data = {header, &records};
		bam_fetch(bam, bam_index, tid, 7000, 7008, &data, fetch_function);
		// sanity check
		expect_true(records.size() == 1);
		std::vector<StartEndPair> blk = get_blocks(records[0]);
		expect_true(blk.size() == 1);
		expect_true(blk[0] == (StartEndPair {7007, 7083}));


		// a check with 409 reads
		std::vector<BAMRecord> recs;
		DataStruct data2 = {header, &recs};
		bam_fetch(bam, bam_index, tid, 7000, 7025, &data2, fetch_function);
		bam_close(bam);	

		// the real values of blocks
		std::vector<std::vector<StartEndPair>> res_blocks {
			{ StartEndPair {7007, 7083} }, 
			{ StartEndPair {7009, 7085} }, 
			{ StartEndPair {7009, 7085} }, 
			{ StartEndPair {7010, 7086} }, 
			{ StartEndPair {7010, 7086} }, 
			{ StartEndPair {7010, 7086} }, 
			{ StartEndPair {7011, 7087} }, 
			{ StartEndPair {7011, 7087} }, 
			{ StartEndPair {7012, 7088} }, 
			{ StartEndPair {7012, 7088} }, 
			{ StartEndPair {7012, 7088} }, 
			{ StartEndPair {7014, 7090} }, 
			{ StartEndPair {7014, 7090} }, 
			{ StartEndPair {7014, 7090} }, 
			{ StartEndPair {7015, 7091} }, 
			{ StartEndPair {7015, 7091} }, 
			{ StartEndPair {7015, 7091} }, 
			{ StartEndPair {7015, 7091} }, 
			{ StartEndPair {7016, 7092} }, 
			{ StartEndPair {7016, 7092} }, 
			{ StartEndPair {7016, 7092} }, 
			{ StartEndPair {7016, 7092} }, 
			{ StartEndPair {7016, 7092} }, 
			{ StartEndPair {7016, 7092} }, 
			{ StartEndPair {7016, 7092} }, 
			{ StartEndPair {7016, 7092} }, 
			{ StartEndPair {7016, 7092} }, 
			{ StartEndPair {7017, 7093} }, 
			{ StartEndPair {7017, 7093} }, 
			{ StartEndPair {7017, 7093} }, 
			{ StartEndPair {7017, 7093} }, 
			{ StartEndPair {7017, 7093} }, 
			{ StartEndPair {7017, 7093} }, 
			{ StartEndPair {7017, 7093} }, 
			{ StartEndPair {7017, 7093} }, 
			{ StartEndPair {7017, 7093} }, 
			{ StartEndPair {7018, 7094} }, 
			{ StartEndPair {7018, 7094} }, 
			{ StartEndPair {7018, 7094} }, 
			{ StartEndPair {7018, 7094} }, 
			{ StartEndPair {7018, 7094} }, 
			{ StartEndPair {7018, 7094} }, 
			{ StartEndPair {7018, 7094} }, 
			{ StartEndPair {7018, 7094} }, 
			{ StartEndPair {7019, 7095} }, 
			{ StartEndPair {7019, 7095} }, 
			{ StartEndPair {7019, 7095} }, 
			{ StartEndPair {7019, 7095} }, 
			{ StartEndPair {7019, 7095} }, 
			{ StartEndPair {7019, 7095} }, 
			{ StartEndPair {7019, 7095} }, 
			{ StartEndPair {7019, 7095} }, 
			{ StartEndPair {7019, 7095} }, 
			{ StartEndPair {7019, 7095} }, 
			{ StartEndPair {7019, 7095} }, 
			{ StartEndPair {7019, 7095} }, 
			{ StartEndPair {7019, 7095} }, 
			{ StartEndPair {7019, 7095} }, 
			{ StartEndPair {7019, 7095} }, 
			{ StartEndPair {7019, 7095} }, 
			{ StartEndPair {7019, 7095} }, 
			{ StartEndPair {7019, 7095} }, 
			{ StartEndPair {7019, 7095} }, 
			{ StartEndPair {7019, 7095} }, 
			{ StartEndPair {7019, 7095} }, 
			{ StartEndPair {7019, 7095} }, 
			{ StartEndPair {7019, 7095} }, 
			{ StartEndPair {7019, 7095} }, 
			{ StartEndPair {7019, 7095} }, 
			{ StartEndPair {7019, 7095} }, 
			{ StartEndPair {7019, 7095} }, 
			{ StartEndPair {7019, 7095} }, 
			{ StartEndPair {7019, 7095} }, 
			{ StartEndPair {7019, 7095} }, 
			{ StartEndPair {7019, 7095} }, 
			{ StartEndPair {7019, 7095} }, 
			{ StartEndPair {7019, 7095} }, 
			{ StartEndPair {7019, 7095} }, 
			{ StartEndPair {7019, 7095} }, 
			{ StartEndPair {7019, 7095} }, 
			{ StartEndPair {7019, 7095} }, 
			{ StartEndPair {7019, 7095} }, 
			{ StartEndPair {7019, 7095} }, 
			{ StartEndPair {7019, 7095} }, 
			{ StartEndPair {7019, 7095} }, 
			{ StartEndPair {7019, 7095} }, 
			{ StartEndPair {7020, 7096} }, 
			{ StartEndPair {7020, 7096} }, 
			{ StartEndPair {7020, 7096} }, 
			{ StartEndPair {7020, 7096} }, 
			{ StartEndPair {7020, 7096} }, 
			{ StartEndPair {7020, 7096} }, 
			{ StartEndPair {7020, 7096} }, 
			{ StartEndPair {7020, 7096} }, 
			{ StartEndPair {7020, 7096} }, 
			{ StartEndPair {7020, 7096} }, 
			{ StartEndPair {7020, 7096} }, 
			{ StartEndPair {7020, 7096} }, 
			{ StartEndPair {7020, 7096} }, 
			{ StartEndPair {7020, 7096} }, 
			{ StartEndPair {7020, 7096} }, 
			{ StartEndPair {7020, 7096} }, 
			{ StartEndPair {7020, 7096} }, 
			{ StartEndPair {7020, 7096} }, 
			{ StartEndPair {7020, 7096} }, 
			{ StartEndPair {7020, 7096} }, 
			{ StartEndPair {7020, 7096} }, 
			{ StartEndPair {7020, 7096} }, 
			{ StartEndPair {7020, 7096} }, 
			{ StartEndPair {7020, 7096} }, 
			{ StartEndPair {7020, 7096} }, 
			{ StartEndPair {7020, 7096} }, 
			{ StartEndPair {7020, 7096} }, 
			{ StartEndPair {7020, 7096} }, 
			{ StartEndPair {7020, 7096} }, 
			{ StartEndPair {7020, 7096} }, 
			{ StartEndPair {7020, 7096} }, 
			{ StartEndPair {7020, 7096} }, 
			{ StartEndPair {7020, 7096} }, 
			{ StartEndPair {7020, 7096} }, 
			{ StartEndPair {7020, 7096} }, 
			{ StartEndPair {7020, 7096} }, 
			{ StartEndPair {7020, 7096} }, 
			{ StartEndPair {7020, 7096} }, 
			{ StartEndPair {7020, 7096} }, 
			{ StartEndPair {7020, 7096} }, 
			{ StartEndPair {7020, 7096} }, 
			{ StartEndPair {7020, 7096} }, 
			{ StartEndPair {7020, 7096} }, 
			{ StartEndPair {7020, 7096} }, 
			{ StartEndPair {7020, 7096} }, 
			{ StartEndPair {7020, 7096} }, 
			{ StartEndPair {7020, 7096} }, 
			{ StartEndPair {7020, 7096} }, 
			{ StartEndPair {7020, 7096} }, 
			{ StartEndPair {7020, 7096} }, 
			{ StartEndPair {7020, 7096} }, 
			{ StartEndPair {7020, 7096} }, 
			{ StartEndPair {7020, 7096} }, 
			{ StartEndPair {7020, 7096} }, 
			{ StartEndPair {7020, 7096} }, 
			{ StartEndPair {7020, 7096} }, 
			{ StartEndPair {7020, 7096} }, 
			{ StartEndPair {7020, 7096} }, 
			{ StartEndPair {7020, 7096} }, 
			{ StartEndPair {7020, 7096} }, 
			{ StartEndPair {7020, 7096} }, 
			{ StartEndPair {7020, 7096} }, 
			{ StartEndPair {7020, 7096} }, 
			{ StartEndPair {7020, 7096} }, 
			{ StartEndPair {7020, 7096} }, 
			{ StartEndPair {7020, 7096} }, 
			{ StartEndPair {7020, 7096} }, 
			{ StartEndPair {7020, 7096} }, 
			{ StartEndPair {7020, 7096} }, 
			{ StartEndPair {7020, 7096} }, 
			{ StartEndPair {7020, 7096} }, 
			{ StartEndPair {7020, 7096} }, 
			{ StartEndPair {7020, 7096} }, 
			{ StartEndPair {7020, 7096} }, 
			{ StartEndPair {7020, 7096} }, 
			{ StartEndPair {7020, 7096} }, 
			{ StartEndPair {7020, 7096} }, 
			{ StartEndPair {7021, 7097} }, 
			{ StartEndPair {7021, 7097} }, 
			{ StartEndPair {7021, 7097} }, 
			{ StartEndPair {7021, 7097} }, 
			{ StartEndPair {7021, 7097} }, 
			{ StartEndPair {7021, 7097} }, 
			{ StartEndPair {7021, 7097} }, 
			{ StartEndPair {7021, 7097} }, 
			{ StartEndPair {7021, 7097} }, 
			{ StartEndPair {7022, 7098} }, 
			{ StartEndPair {7022, 7098} }, 
			{ StartEndPair {7022, 7098} }, 
			{ StartEndPair {7022, 7098} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7077}, StartEndPair {7077, 7097} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7023, 7099} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7099} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7099} }, 
			{ StartEndPair {7024, 7099} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7099} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7099} }, 
			{ StartEndPair {7024, 7099} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7099} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7099} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7096} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7099} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7099} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7098} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7099} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7098} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7099} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7099} }, 
			{ StartEndPair {7024, 7099} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7099} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7099} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7100} }, 
			{ StartEndPair {7024, 7100} }
		};

		std::vector<std::vector<StartEndPair>> rec_blocks = map<BAMRecord, std::vector<StartEndPair>>(recs, [](const BAMRecord &b) { return get_blocks(b); });
	
		expect_true(res_blocks.size() == rec_blocks.size());
		for (int i = 0; i < (int)res_blocks.size(); i++) {
			expect_true(compare_stream(res_blocks[i], rec_blocks[i]));
		}
	}

	test_that("full group_bam2isoform function executes correctly") {
// from sc_longread import *
// from parse_gene_anno import *
// from parse_config import parse_json_config
// bam = '/Users/voogd.o/Documents/FLAMESData/data/align2genome.bam'
// genomefa = '/Volumes/voogd.o/FLAMES/inst/extdata/SIRV_genomefa.fasta'
// gff3 = '/Volumes/voogd.o/FLAMES/inst/extdata/SIRV_anno.gtf'
// isoform_gff3 = '/Users/voogd.o/Documents/FLAMESintermediate/isoform_annotated.gff3'
// config = '/Volumes/voogd.o/FLAMES/inst/extdata/SIRV_config_default.json'
// tss_tes_stat = '/Users/voogd.o/Documents/FLAMESintermediate/tss_tes.bedgraph'
// chr_to_gene, transcript_dict, gene_to_transcript, transcript_to_exon = parse_gff_tree(gff3)
// transcript_to_junctions = {tr: blocks_to_junctions(transcript_to_exon[tr]) for tr in transcript_to_exon}
// remove_similar_tr(transcript_dict, gene_to_transcript, transcript_to_exon)
// gene_dict = get_gene_flat(gene_to_transcript, transcript_to_exon)
// chr_to_blocks = get_gene_blocks(gene_dict, chr_to_gene, gene_to_transcript)
// config_dict = parse_json_config(config)
// group_bam2isoform(bam, isoform_gff3, tss_tes_stat, "", chr_to_blocks, gene_dict, transcript_to_junctions, transcript_dict, genomefa, config=config_dict["isoform_parameters"], downsample_ratio=1,raw_gff3=None)

		// Function downloadFile("download.file");
		// std::string bam_in = downloadFile(_["url"]="https://raw.githubusercontent.com/OliverVoogd/FLAMESData/master/data/align2genome.bam", _["destfile"]=get_tempfile(".bam"));
		// std::string out_gff3 = get_tempfile("isoform_annotated.gff3"); 
		// std::string out_stat = get_tempfile("tss_tes.bedgraph");
		// std::string fa_f = get_extdata("SIRV_genomefa.fasta");
		// IsoformParameters isoform_parameters;
		// std::string raw_gff3;

		// std::string gff3 = get_extdata("SIRV_anno.gtf");
		// GeneAnnoParser parser (gff3);
		// GFFData data = parser.parse();

		// // data manipulation
		// std::unordered_map<std::string, Junctions> transcript_to_junctions = 
		// 	map<std::string, std::vector<StartEndPair>, Junctions>(data.transcript_to_exon, blocks_to_junctions);

		// remove_similar_tr(data.gene_to_transcript, data.transcript_to_exon, 10);
		// std::unordered_map<std::string, std::vector<StartEndPair>> gene_dict =
		// 	get_gene_flat(data.gene_to_transcript, data.transcript_to_exon);
		// std::unordered_map<std::string, std::vector<GeneBlocks>> chr_to_blocks = 
		// 	get_gene_blocks(gene_dict, data.chr_to_gene, data.gene_to_transcript);
		// group_bam2isoform(
		// 	bam_in, out_gff3, out_stat, 
		// 	chr_to_blocks, gene_dict, 
		// 	transcript_to_junctions, data.gene_to_transcript,
		// 	fa_f, isoform_parameters, raw_gff3);
	}
}

// [[Rcpp::export]]
void what2() {
	Rcpp::Function downloadFile("download.file");
	std::string bam_dir = get_tempdir();
	std::string bam_in = bam_dir + std::string("/align2genome.bam");
	std::string bam_idx_in = bam_dir + std::string("/align2genome.bam.bai");
	int bam_in_flag = Rcpp::as<int> (downloadFile(Rcpp::_["url"]="https://raw.githubusercontent.com/OliverVoogd/FLAMESData/master/data/align2genome.bam", Rcpp::_["destfile"]=bam_in));
	int bam_idx_flag = Rcpp::as<int> (downloadFile(Rcpp::_["url"]="https://raw.githubusercontent.com/OliverVoogd/FLAMESData/master/data/align2genome.bam.bai", Rcpp::_["destfile"]=bam_idx_in));
	if (bam_in_flag != 0 || bam_idx_flag != 0) {
		Rcpp::Rcout << "downloaded failed";
		return;
	}

	std::string out_gff3 = get_tempfile("isoform_annotated.gff3"); 
	std::string out_stat = get_tempfile("tss_tes.bedgraph");
	Rcpp::Rcout << "Out GFF3: " << out_gff3 << "\n";
	Rcpp::Rcout << "Out STAT: " << out_stat << "\n";
	std::string fa_f = get_extdata("SIRV_genomefa.fasta");
	IsoformParameters isoform_parameters;
	std::string raw_gff3;

	std::string gff3 = get_extdata("SIRV_anno.gtf");
	GeneAnnoParser parser (gff3);
	GFFData data = parser.parse();

	// data manipulation
	std::unordered_map<std::string, Junctions> transcript_to_junctions = 
		map<std::string, std::vector<StartEndPair>, Junctions>(data.transcript_to_exon, blocks_to_junctions);

	remove_similar_tr(data.gene_to_transcript, data.transcript_to_exon, 10);
	std::unordered_map<std::string, std::vector<StartEndPair>> gene_dict =
		get_gene_flat(data.gene_to_transcript, data.transcript_to_exon);
	std::unordered_map<std::string, std::vector<GeneBlocks>> chr_to_blocks = 
		get_gene_blocks(gene_dict, data.chr_to_gene, data.gene_to_transcript);
	
	group_bam2isoform(
		bam_in, out_gff3, out_stat, 
		chr_to_blocks, gene_dict, 
		transcript_to_junctions, data.transcript_dict,
		fa_f, isoform_parameters, raw_gff3);
}	
