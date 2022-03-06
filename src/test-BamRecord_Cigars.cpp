#include <string>
#include <vector>
#include <sstream>
#include <algorithm>

#include <testthat.h>

#include "classes/BamRecord.h"
#include "classes/DataStruct.h"
#include "utility/cigars.h"
#include "utility/bam.h"
#include "test_utilities.h"

class CigarMap {
private:
	std::unordered_map<std::string, std::vector<CigarPair>> cigars;
public:
	void add(const char *qname, const std::vector<CigarPair> &cig) {
		if (!contains(qname)) {
			cigars[std::string(qname)] = cig;
		}
	}

	bool contains(const char *qname) {
		return cigars.count(std::string(qname)) > 0;
	}

	int size() {
		return cigars.size();
	}

	std::unordered_map<std::string, std::vector<CigarPair>> get() {
		return cigars;
	}

	std::vector<CigarPair> at(const char *qname) {
		if (contains(qname)) {
			return cigars.at(std::string(qname));
		} else {
			return std::vector<CigarPair>();
		}
	}

	CigarMap() {
		cigars = std::unordered_map<std::string, std::vector<CigarPair>>();
	}
};

struct FetchStruct {
	const bam_header_t *header;
	BAMRecord *rec;
};

// static inline std::string printCigarPair(const CigarPair &cp) {
// 	return "(" + std::to_string(cp.op) + std::string(",") + std::to_string(cp.len) + std::string(")");
// }

// static inline std::string printCigar(std::vector<CigarPair> m) {
// 	return "[" 
// 		+ std::accumulate(m.begin()+1, m.end(), printCigarPair(*m.begin()), 
// 			[](const std::string &acc, const CigarPair &cur){ return acc + "," + printCigarPair(cur); })
// 		+ "]";
// }

static inline bool compareBAMRecord(const BAMRecord &a, const BAMRecord &b) {
	return 
		a.cigar_string == b.cigar_string &&
		a.reference_name == b.reference_name && 
		a.reference_start == b.reference_start && 
		a.reference_end == b.reference_end && 
		a.read_name == b.read_name && 
		a.AS_tag == b.AS_tag && 
		a.query_alignment_length == b.query_alignment_length && 
		a.mapping_quality == b.mapping_quality; // && 
		// a.flag == b.flag;
	// should we compare cigar vector?
}

// fetch function for testing generate_cigar_pairs over a number of entries
static int testFetch1(const bam1_t *b, void *data) {
	CigarMap *cigar = (CigarMap *)data;

	const char *qname = bam1_qname(b);
	std::vector<CigarPair> pairs = generate_cigar_pairs(b);
	cigar->add(qname, pairs);
	
	return 1;
}

// fetch function for a single BAMRecord entry
static int testFetch2(const bam1_t *b, void *data) {
	FetchStruct *fs = (FetchStruct *)data;

	const bam_header_t *header = fs->header;
	BAMRecord *rec = fs->rec;

	if (std::string(bam1_qname(b)) == "ERR032227.10543296") {
		*rec = read_record(b, header);
	}
	return 1;
}

// fetch function for testing a large number of BAMRecord entries
static int fetch_function(const bam1_t *b, void *data)
{
    DataStruct * data_struct = (DataStruct*)data;

    std::vector<BAMRecord> * records = data_struct->records;
    bam_header_t * header = data_struct->header;

    BAMRecord rec = read_record(b, header);
    records->push_back(rec);
	return 0;
}

context("BamRecord & Utility BAM/CIGAR Function Testing") {
	test_that("Generated cigar pairs are the same as pysam") {
		// bam = pysam.AlignmentFile("/Users/voogd.o/Downloads/Test1-ready.bam");
		// reg = [x for x in bam.fetch("chrM", 7000, 7008)][0]
		std::string sBamIn = get_extdata("Test1-ready.bam");
		const char *bam_in = sBamIn.c_str();
		bamFile bam = bam_open(bam_in, "r");
		bam_index_t *bam_index = bam_index_load(bam_in);
		int tid = 0;

		CigarMap cigar;
		bam_fetch(bam, bam_index, tid, 7000, 7008, &cigar, testFetch1);
		

		// the resulting cigar vec should be {(0, 76)}
		expect_true(cigar.contains("ERR032227.10543296"));
		if (cigar.contains("ERR032227.10543296")) {
			expect_true(cigar.at("ERR032227.10543296").size() == 1);
			std::vector<CigarPair> res {CigarPair{0, 76}};
			expect_true(compare_unordered(cigar.at("ERR032227.10543296"), res));
		}

		// compare another position
		CigarMap cigar2;
		bam_fetch(bam, bam_index, tid, 7050, 7060, &cigar2, testFetch1);
		bam_close(bam);
		
		expect_true(cigar2.contains("ERR032227.11196628"));
		if (cigar2.contains("ERR032227.11196628")) {
			expect_true(cigar2.at("ERR032227.11196628").size() == 3);
			std::vector<CigarPair> res {CigarPair{0, 54}, CigarPair{1,1}, CigarPair{0, 20}};
			expect_true(compare_unordered(cigar2.at("ERR032227.11196628"), res));
		}
	}

	test_that("generated cigar strings are correct") {
		std::string sBamIn = get_extdata("Test1-ready.bam");
		const char *bam_in = sBamIn.c_str();
		bamFile bam = bam_open(bam_in, "r");
		bam_index_t *bam_index = bam_index_load(bam_in);
		int tid = 0;

		CigarMap cigar;
		bam_fetch(bam, bam_index, tid, 7000, 7008, &cigar, testFetch1);

		expect_true(generate_cigar(cigar.at("ERR032227.10543296")) == "76M");

		CigarMap cigar2;
		bam_fetch(bam, bam_index, tid, 7050, 7060, &cigar2, testFetch1);
		expect_true(generate_cigar(cigar2.at("ERR032227.11196628")) == "54M1I20M");
		expect_true(generate_cigar(cigar2.at("ERR032227.11204441")) == "54M1D22M");
		bam_close(bam);
	}

	test_that("cigars can be smoothed") {
			std::vector<CigarPair> cigar {
				{0, 7}, {2, 1}, {0, 4}, {2, 1}, {0, 27}, {1, 1}, {0, 9}, {2, 1}, {0, 4}, {2, 1}, {0, 3}, {2, 2}, {0, 3}, {2, 1}, {0, 17}, {1, 1}, {0, 2}, {1, 2}, {0, 28}, {2, 1}, {0, 14}, {2, 1}, {0, 8}, {3, 838}, {0, 15}, {2, 1}, {0, 11}, {2, 1}, {0, 18}, {3, 86}, {0, 19}, {2, 1}, {0, 15}, {1, 1}, {0, 2}, {3, 114}, {0, 17}, {2, 1}, {0, 7}, {2, 1}, {0, 7}, {1, 1}, {0, 3}, {2, 1}, {0, 8}, {3, 983}, {0, 15}, {2, 5}, {0, 7}, {2, 1}, {0, 60}, {2, 3}, {0, 15}, {3, 79}, {0, 12}, {2, 1}, {0, 17}, {2, 1}, {0, 5}, {2, 1}, {0, 13}, {1, 1}, {0, 11}, {2, 1}, {0, 37}, {2, 1}, {0, 6}, {2, 3}, {0, 1}, {2, 1}, {0, 2}, {2, 1}, {0, 3}, {2, 1}, {0, 13}, {2, 3}, {0, 6}, {2, 1}, {0, 2}, {2, 2}, {0, 3}, {1, 1}, {0, 12}, {3, 1737}, {0, 9}, {1, 1}, {0, 10}, {2, 1}, {0, 7}, {2, 1}, {0, 2}, {1, 2}, {0, 9}, {2, 1}, {0, 20}, {1, 1}, {0, 6}, {1, 1}, {0, 4}, {3, 93}, {0, 14}, {1, 5}, {0, 5}, {1, 3}, {0, 9}, {2, 1}, {0, 4}, {1, 2}, {0, 11}, {2, 1}, {0, 15}, {2, 1}, {0, 14}, {2, 1}, {0, 1}, {2, 2}, {0, 4}, {3, 483}, {0, 7}, {1, 1}, {0, 12}, {2, 1}, {0, 34}, {2, 2}, {0, 4}, {3, 158}, {0, 2}, {1, 1}, {0, 7}, {1, 1}, {0, 3}, {2, 2}, {0, 10}, {2, 1}, {0, 19}, {1, 1}, {0, 3}, {1, 2}, {0, 12}, {2, 4}, {0, 5}, {2, 1}, {0, 7}, {2, 1}, {0, 9}, {2, 4}, {0, 3}, {2, 1}, {0, 1}, {2, 2}, {0, 2}, {1, 1}, {0, 10}, {1, 1}, {0, 5}, {2, 1}, {0, 5}, {2, 2}, {0, 3}, {3, 374}, {0, 15}, {2, 3}, {0, 6}, {2, 3}, {0, 22}, {1, 1}, {0, 9}, {1, 1}, {0, 4}, {2, 4}, {0, 4}, {2, 2}, {0, 16}, {2, 1}, {0, 2}, {2, 1}, {0, 7}, {1, 2}, {0, 1}, {2, 1}, {0, 13}, {1, 2}, {0, 16}, {2, 1}, {3, 187}, {0, 9}, {1, 2}, {0, 2}, {2, 1}, {0, 3}, {2, 2}, {0, 16}, {1, 1}, {0, 3}, {2, 2}, {0, 12}, {2, 3}, {0, 11}, {2, 1}, {0, 6}, {1, 2}, {0, 1}, {1, 1}, {0, 4}, {1, 4}, {0, 14}, {2, 1}, {0, 14}, {2, 1}, {0, 2}, {1, 3}, {0, 9}, {2, 1}, {0, 16}, {2, 1}, {0, 28}, {3, 374}, {0, 36}, {2, 3}, {0, 1}, {2, 1}, {0, 2}, {2, 2}, {0, 21}, {1, 1}, {0, 3}, {2, 1}, {0, 2}, {1, 3}, {0, 4}, {2, 1}, {0, 4}, {3, 136}, {0, 7}, {2, 1}, {0, 1}, {2, 2}, {0, 5}, {2, 1}, {0, 8}, {2, 1}, {0, 6}, {2, 1}, {0, 5}, {2, 2}, {0, 15}, {2, 3}, {0, 8}, {2, 1}, {0, 18}, {2, 1}, {0, 11}, {2, 2}, {0, 7}, {2, 7}, {0, 1}, {2, 1}, {0, 8}, {2, 1}, {0, 7}, {2, 6}, {0, 13}, {2, 4}, {0, 9}, {2, 7}, {0, 1}, {2, 1}, {0, 4}, {2, 4}, {0, 11}, {1, 1}, {0, 8}, {2, 1}, {0, 20}, {2, 4}, {0, 12}, {2, 1}, {0, 5}, {2, 3}, {0, 18}, {2, 1}, {0, 2}, {2, 2}, {0, 5}, {2, 3}, {0, 20}, {2, 2}, {0, 6}, {2, 2}, {0, 4}, {1, 1}, {0, 9}, {2, 2}, {0, 12}, {2, 4}, {0, 4}, {2, 1}, {0, 16}, {2, 2}, {0, 10}, {2, 2}, {0, 2}, {2, 1}, {0, 1}, {2, 1}, {0, 9}, {2, 2}, {0, 20}, {1, 1}, {0, 22}, {2, 1}, {0, 3}, {2, 7}, {0, 25}, {2, 5}, {0, 12}, {3, 73}, {0, 13}, {2, 1}, {0, 6}, {2, 2}, {0, 20}, {2, 1}, {0, 1}, {2, 1}, {0, 3}, {2, 1}, {0, 22}, {2, 2}, {0, 10}, {2, 1}, {0, 37}, {2, 3}, {0, 7}, {3, 2273}, {0, 24}, {2, 2}, {0, 14}, {2, 5}, {0, 10}, {2, 1}, {0, 9}, {1, 1}, {0, 2}, {2, 1}, {0, 3}, {2, 1}, {0, 16}, {1, 1}, {0, 6}, {2, 3}, {0, 5}, {2, 5}, {0, 13}, {2, 1}, {0, 10}, {4, 55}
			};

			std::vector<CigarPair> smoothed {
				{0, 135}, {3, 838}, {0, 46}, {3, 86}, {0, 37}, {3, 114}, {0, 45}, {3, 983}, {0, 106}, {3, 79}, {0, 160}, {3, 1737}, {0, 70}, {3, 93}, {0, 83}, {3, 483}, {0, 60}, {3, 158}, {0, 125}, {3, 374}, {0, 131}, {3, 187}, {0, 163}, {3, 374}, {0, 81}, {3, 136}, {0, 483}, {3, 73}, {0, 131}, {3, 2273}, {0, 131}, {4, 55}
			};

			std::vector<CigarPair> res = smooth_cigar(cigar, 10);

			expect_true(smoothed.size() == res.size());
			expect_true(compare_stream(smoothed, res));
	}

	test_that("BAMRecord is generated correctly from bam1_t") {
		std::string sBamIn = get_extdata("Test1-ready.bam");
		const char *bam_in = sBamIn.c_str();
		bamFile bam = bam_open(bam_in, "r");
		bam_index_t *bam_index = bam_index_load(bam_in);
		bam_header_t *header = bam_header_read(bam);
		int tid = 0;

		BAMRecord rec;
		FetchStruct fs {header, &rec};

		bam_fetch(bam, bam_index, tid, 7000, 7008, &fs, testFetch2);		
		BAMRecord res { std::vector<CigarPair>{ {0, 76} }, "76M", "chrM", 7007, 7083, "ERR032227.10543296", 0, 76, 50, read_flag(99)};
		
		expect_true(compareBAMRecord(rec, res));


		bam_close(bam);
	}

	test_that("Check correct reading of bam for large number of reads") {
		std::string bamIn = get_extdata("Test1-ready.bam");
		bamFile bam = bam_open(bamIn.c_str(), "r");
		bam_index_t *bam_index = bam_index_load(bamIn.c_str());
		bam_header_t *header = bam_header_read(bam);
		int tid = bam_get_tid(header, "chrM");
		// do a check on a position with more records
		std::vector<BAMRecord> recs;
		DataStruct data2 = {header, &recs};
		bam_fetch(bam, bam_index, tid, 7000, 7025, &data2, fetch_function);
		bam_close(bam);	

		// recs now contains 409 bam records which we can check against a few key values
		std::vector<std::string> res_qnames {
			"ERR032227.10543296", "ERR032227.10212620", "ERR032227.10782111", "ERR032227.10060210", "ERR032227.11020015", "ERR032227.11305478", "ERR032227.10748414", "ERR032227.11283115", "ERR032227.10311580", "ERR032227.10533999", "ERR032227.11278070", "ERR032227.1009846", "ERR032227.10113962", "ERR032227.10898858", "ERR032227.1004669", "ERR032227.1048664", "ERR032227.10548216", "ERR032227.10930732", "ERR032227.10023673", "ERR032227.10042289", "ERR032227.10193728", "ERR032227.10532659", "ERR032227.10707595", "ERR032227.10757488", "ERR032227.10909480", "ERR032227.10954055", "ERR032227.1102213", "ERR032227.10320298", "ERR032227.10351125", "ERR032227.10501899", "ERR032227.10589856", "ERR032227.10785442", "ERR032227.10805014", "ERR032227.11120114", "ERR032227.11348276", "ERR032227.11358297", "ERR032227.10167843", "ERR032227.10178286", 
			"ERR032227.10403142", "ERR032227.10483558", "ERR032227.10514960", "ERR032227.10540949", "ERR032227.10907484", "ERR032227.1139107", "ERR032227.1004880", "ERR032227.10128936", "ERR032227.10164143", "ERR032227.10187615", "ERR032227.10220095", "ERR032227.10225741", "ERR032227.10274819", "ERR032227.10322564", "ERR032227.10341548", "ERR032227.10347198", "ERR032227.103743", "ERR032227.10407341", "ERR032227.10458683", "ERR032227.10486139", "ERR032227.10497590", "ERR032227.10512576", "ERR032227.10543585", "ERR032227.10545891", "ERR032227.10569341", "ERR032227.10570150", "ERR032227.10603169", "ERR032227.10616675", "ERR032227.10768707", "ERR032227.10805568", "ERR032227.10839198", "ERR032227.10846281", "ERR032227.10944860", "ERR032227.1094685", "ERR032227.10949965", "ERR032227.10995599", "ERR032227.11008542", "ERR032227.11044617", 
			"ERR032227.11053135", "ERR032227.1112306", "ERR032227.1115690", "ERR032227.11190484", "ERR032227.1120433", "ERR032227.11245650", "ERR032227.1131244", "ERR032227.11325606", "ERR032227.11355838", "ERR032227.11391123", "ERR032227.10015612", "ERR032227.1001892", "ERR032227.10067714", "ERR032227.10161514", "ERR032227.10186085", "ERR032227.10194731", "ERR032227.10227045", "ERR032227.10248970", "ERR032227.10251051", "ERR032227.10267256", "ERR032227.10317321", "ERR032227.10327370", "ERR032227.10346923", "ERR032227.10411552", "ERR032227.10417427", "ERR032227.10428470", "ERR032227.10466280", "ERR032227.10468909", "ERR032227.10509552", "ERR032227.10514399", "ERR032227.10521093", "ERR032227.10522925", "ERR032227.1055273", "ERR032227.10555534", "ERR032227.10556798", "ERR032227.10558381", "ERR032227.10628069", "ERR032227.10633261", 
			"ERR032227.10639590", "ERR032227.10641781", "ERR032227.10758186", "ERR032227.10770439", "ERR032227.10779135", "ERR032227.10785127", "ERR032227.10792627", "ERR032227.10835148", "ERR032227.10853738", "ERR032227.10858581", "ERR032227.10871989", "ERR032227.10887884", "ERR032227.10894473", "ERR032227.10903091", "ERR032227.10914756", "ERR032227.10924006", "ERR032227.10947028", "ERR032227.10975703", "ERR032227.11019321", "ERR032227.11020069", "ERR032227.11024278", "ERR032227.11036120", "ERR032227.11053885", "ERR032227.11056753", "ERR032227.11111853", "ERR032227.11130054", "ERR032227.11145871", "ERR032227.11169828", "ERR032227.11181354", "ERR032227.1119300", "ERR032227.11193379", "ERR032227.11203118", "ERR032227.11222431", "ERR032227.11233205", "ERR032227.11270957", "ERR032227.11292359", "ERR032227.11304234", "ERR032227.11305461", 
			"ERR032227.11308208", "ERR032227.11318467", "ERR032227.11321474", "ERR032227.11327104", "ERR032227.1132846", "ERR032227.11329094", "ERR032227.1135841", "ERR032227.11365968", "ERR032227.11376940", "ERR032227.11403060", "ERR032227.11415682", "ERR032227.10340164", "ERR032227.10605326", "ERR032227.10810673", "ERR032227.10823886", "ERR032227.11020294", "ERR032227.11058936", "ERR032227.11068243", "ERR032227.11159103", "ERR032227.11221924", "ERR032227.1000832", "ERR032227.10715782", "ERR032227.11218866", "ERR032227.11258017", "ERR032227.10005617", "ERR032227.10040273", "ERR032227.10040978", "ERR032227.10049464", "ERR032227.10053208", "ERR032227.10056048", "ERR032227.10062266", "ERR032227.1007294", "ERR032227.1009594", "ERR032227.10117932", "ERR032227.10121955", "ERR032227.10133874", "ERR032227.10139264", "ERR032227.1013953", "ERR032227.10140593", 
			"ERR032227.10156979", "ERR032227.10207697", "ERR032227.10224254", "ERR032227.10228976", "ERR032227.10233155", "ERR032227.10276459", "ERR032227.10278940", "ERR032227.10280405", "ERR032227.10295438", "ERR032227.10304513", "ERR032227.10312061", "ERR032227.10337006", "ERR032227.10365556", "ERR032227.10381933", "ERR032227.10383330", "ERR032227.10408984", "ERR032227.10421403", "ERR032227.10449470", "ERR032227.10451096", "ERR032227.10491068", "ERR032227.10497035", "ERR032227.1051201", "ERR032227.10532276", "ERR032227.10535379", "ERR032227.10544932", "ERR032227.10565160", "ERR032227.10575788", "ERR032227.10591709", "ERR032227.10600607", "ERR032227.10622313", "ERR032227.10645676", "ERR032227.10666870", "ERR032227.10667897", "ERR032227.10668557", "ERR032227.10672382", "ERR032227.10678812", "ERR032227.10688359", "ERR032227.1070109", "ERR032227.107176", 
			"ERR032227.10720185", "ERR032227.10730221", "ERR032227.10738971", "ERR032227.10760190", "ERR032227.10814454", "ERR032227.10827569", "ERR032227.10841862", "ERR032227.10843256", "ERR032227.10888375", "ERR032227.10890793", "ERR032227.1089282", "ERR032227.10897661", "ERR032227.10905853", "ERR032227.10906185", "ERR032227.10930899", "ERR032227.10941393", "ERR032227.10942498", "ERR032227.10943145", "ERR032227.1096683", "ERR032227.10986576", "ERR032227.10992072", "ERR032227.11001288", "ERR032227.1100576", "ERR032227.11008040", "ERR032227.11032848", "ERR032227.1106370", "ERR032227.11074293", "ERR032227.11080311", "ERR032227.11084710", "ERR032227.11118450", "ERR032227.11123474", "ERR032227.11127581", "ERR032227.11173856", "ERR032227.11196628", "ERR032227.11197117", "ERR032227.11199439", "ERR032227.11199901", "ERR032227.11226522", "ERR032227.11239404", 
			"ERR032227.11241721", "ERR032227.11250328", "ERR032227.11258764", "ERR032227.11259311", "ERR032227.11270515", "ERR032227.11298361", "ERR032227.11304346", "ERR032227.11332661", "ERR032227.11333780", "ERR032227.11379556", "ERR032227.11400998", "ERR032227.11404948", "ERR032227.10001146", "ERR032227.10013108", "ERR032227.1003007", "ERR032227.10037703", "ERR032227.1004312", "ERR032227.10049241", "ERR032227.10054692", "ERR032227.10060224", "ERR032227.10072447", "ERR032227.10096129", "ERR032227.10101354", "ERR032227.10104736", "ERR032227.10114987", "ERR032227.10122971", "ERR032227.10132315", "ERR032227.10137331", "ERR032227.10140081", "ERR032227.10140210", "ERR032227.10157999", "ERR032227.10173338", "ERR032227.10182712", "ERR032227.10201823", "ERR032227.10219493", "ERR032227.1022397", "ERR032227.10233869", "ERR032227.10239568", "ERR032227.10254554", 
			"ERR032227.1026788", "ERR032227.10277751", "ERR032227.10285308", "ERR032227.10285446", "ERR032227.10296579", "ERR032227.10298182", "ERR032227.10310683", "ERR032227.10312808", "ERR032227.10332335", "ERR032227.10333468", "ERR032227.10346320", "ERR032227.10347242", "ERR032227.10347682", "ERR032227.10375112", "ERR032227.10387561", "ERR032227.10400862", "ERR032227.10432595", "ERR032227.10454461", "ERR032227.10454784", "ERR032227.10458903", "ERR032227.10462736", "ERR032227.10465705", "ERR032227.10469402", "ERR032227.10526120", "ERR032227.10537658", "ERR032227.10551650", "ERR032227.10555726", "ERR032227.10558587", "ERR032227.10562444", "ERR032227.10568476", "ERR032227.10569663", "ERR032227.10579711", "ERR032227.10609153", "ERR032227.1061030", "ERR032227.10617181", "ERR032227.106317", "ERR032227.10678800", "ERR032227.10681192", "ERR032227.10681233", 
			"ERR032227.10689934", "ERR032227.10695594", "ERR032227.10701751", "ERR032227.10709029", "ERR032227.10719749", "ERR032227.10724745", "ERR032227.1073042", "ERR032227.10776764", "ERR032227.10779161", "ERR032227.10787042", "ERR032227.10787102", "ERR032227.10827760", "ERR032227.10850303", "ERR032227.10867160", "ERR032227.10881711", "ERR032227.10895008", "ERR032227.1095675", "ERR032227.10969064", "ERR032227.10984976", "ERR032227.10985033", "ERR032227.1101175", "ERR032227.11023008", "ERR032227.11025027", "ERR032227.11035070", "ERR032227.11043854", "ERR032227.11053761", "ERR032227.11064531", "ERR032227.11066477", "ERR032227.1107254", "ERR032227.11077880", "ERR032227.11079171", "ERR032227.11099022", "ERR032227.11107295", "ERR032227.11111712", "ERR032227.11112998", "ERR032227.1112849", "ERR032227.11128522", "ERR032227.11133342", "ERR032227.11134514", 
			"ERR032227.11141636", "ERR032227.11162462", "ERR032227.11181577", "ERR032227.11213134", "ERR032227.11227573", "ERR032227.11233162", "ERR032227.11244336", "ERR032227.11258889", "ERR032227.11274880", "ERR032227.11298276", "ERR032227.11302629", "ERR032227.11315409", "ERR032227.11318684", "ERR032227.11324177", "ERR032227.1133195", "ERR032227.11332245", "ERR032227.11353438", "ERR032227.1135960", "ERR032227.113636", "ERR032227.11365355", "ERR032227.11374519", "ERR032227.11414068", "ERR032227.11416277"
		};

		// compare stored read_names (qnames) against to ensure we are dealing with the correct data
		std::vector<std::string> rec_qnames = map<BAMRecord, std::string>(recs, [](const BAMRecord &a) { return a.read_name; });
		expect_true(compare_stream(rec_qnames, res_qnames));
	}
}
