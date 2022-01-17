#include <string>
#include <sstream>
#include <fstream>
#include <filesystem>
#include <map>
#include <set>
#include <vector>
#include <algorithm>
#include <functional>
#include <numeric>
#include <Rcpp.h>

#include "gtf_to_bed.h"
#include "Parser.h"

// inline functions to convert a vector to a comma-separated strings
static inline std::string vector_to_str(std::vector<int> vector) {
	return std::accumulate(
					std::next(vector.begin()), vector.end(), 
					std::to_string(vector[0]), 
					[](std::string a, int b){
						return a + "," + std::to_string(b);
					}
	);
}

// static inline std::string vector_to_str(std::vector<std::string> vector) {
// 	return std::accumulate(
// 					std::next(vector.begin()), vector.end(), 
// 					vector[0], 
// 					[](std::string a, std::string b){
// 						return a + "," + b;
// 					}
// 	);
// }

static inline std::string set_to_str(std::set<std::string> set) {
	return std::accumulate(
					std::next(set.begin()), set.end(),
					*set.begin(),
					[](std::string a, std::string b){
						return a + "," + b;
					}
	);
}

// lambda function to write rows to a bed file
// format described here: https://en.wikipedia.org/wiki/BED_(file_format)
static inline void bed_write_row (
		std::ofstream &file, 
		std::string chrom,
		int chromStart,
		int chromEnd,
		std::string name,
		int score,
		std::string strand,
		int thickStart,
		int thickEnd,
		std::string itemRgb,
		int blockCount,
		std::string blockSizes,
		std::string blockStarts) {
	if (file.is_open()) {
		file << chrom << "\t"
			<< chromStart << "\t"
			<< chromEnd << "\t"
			<< name << "\t"
			<< score << "\t"
			<< strand << "\t"
			<< thickStart << "\t"
			<< thickEnd << "\t"
			<< itemRgb << "\t"
			<< blockCount << "\t"
			<< blockSizes << "\t"
			<< blockStarts << "\n";
	}
};

// lambda function to write rows to a psl file
// format described here: http://genome.ucsc.edu/FAQ/FAQformat#format2
static inline void psl_write_row(
		std::ofstream &bed, 
		int matches,
		int misMatches,
		int repMatches,
		int nCount,
		int qNumInsert,
		int qBaseInsert,
		int tNumInsert,
		int tBaseInsert,
		std::string strand,
		std::string qName,
		int qSize,
		int qStart,
		int qEnd,
		std::string tName,
		int tSize,
		int tStart,
		int tEnd,
		int blockCount,
		std::string blockSizes,
		std::string qStarts,
		std::string tStarts) {

	bed << matches << "\t"
		<< misMatches << "\t"
		<< repMatches << "\t"
		<< nCount << "\t"
		<< qNumInsert << "\t"
		<< qBaseInsert << "\t"
		<< tNumInsert << "\t"
		<< qBaseInsert << "\t"
		<< tNumInsert << "\t"
		<< tBaseInsert << "\t"
		<< strand << "\t"
		<< qName << "\t"
		<< qSize << "\t"
		<< qStart << "\t"
		<< qEnd << "\t"
		<< tName << "\t"
		<< tSize << "\t"
		<< tStart << "\t"
		<< tEnd << "\t"
		<< blockCount << "\t"
		<< blockSizes << "\t"
		<< qStarts << "\t"
		<< tStarts << "\n";
};

// generic function to write any amount of data to a tab-separated CSV row
static void csv_write_row (std::ofstream& file, std::vector<std::string> data) {
	for (auto entry = data.begin(); entry != data.end() -1 ; entry++) {
		file << *entry << '\t';
	}
	file << data.back() << '\n';
};


// [[Rcpp::export]]
void
gtf_to_bed_cpp(std::string in_gtf, std::string out_bed, std::string chrom_sizes_file)
{
	std::ifstream gtf (in_gtf);
	std::ofstream bed (out_bed);
	std::ifstream chrom_sizes;

	bool is_bed = (out_bed.substr(out_bed.length() - 3, out_bed.length()) != "psl") ? true : false;

	if (chrom_sizes_file.length()) {
		chrom_sizes.open(chrom_sizes_file); 
	}

	std::map<std::string, int> chrom_to_size;

	if (chrom_sizes.is_open()) {
		chrom_to_size = parsePairsToMap(chrom_sizes);
	}

	std::set<std::string> missing_chroms;

	std::string prev_transcript;
	std::string prev_chrom;
	std::string prev_gene;
	std::string prev_strand;

	int blockcount;
	std::vector<int> blockstarts;
	std::vector<int> blocksizes;
	std::vector<int> qstarts;

	std::string gtf_line;
	std::string this_transcript;
	int tstart;
	int tend;

	std::string chrom;
	std::string strand;
	int start;
	int end;


	while (std::getline(gtf, gtf_line)) {
		if (gtf_line[0] == '#') {
			continue;
		}
		
		// std::stringstream gtf_line_stream (gtf_line);
		std::vector<std::string> values = parseLine(gtf_line, '\t');
		chrom = values[0];
		std::string ty = values[2]; // only used in this while statement
		start = std::stoi(values[3]) - 1;
		end = std::stoi(values[4]);
		strand = values[6];

		if (ty != "exon") {
			continue;
		}
		// just get a substring of the attributes entry
		this_transcript = values[8].substr(values[8].find("transcript_id") + 15, values[8].length());
		this_transcript = parseUntilChar(this_transcript, '"').first;

    	// once all the exons for a transcript are read, write the psl/bed entry
		if (this_transcript != prev_transcript) {
			if (prev_transcript.length()) {
				blockcount = blockstarts.size();
				if (blockcount > 1 && blockstarts[0] > blockstarts[1]) { // we need to reverse the exons
					std::reverse(blocksizes.begin(), blocksizes.end());
					std::reverse(blockstarts.begin(), blockstarts.end());
				}
				// target (eg chrom)
				tstart = blockstarts.front();
				tend = blockstarts.back() + blocksizes.back();

				// query (eg transcript)
				int qsize = std::accumulate(blocksizes.begin(), blocksizes.end(), 0);

				std::string qname = prev_transcript; 
				
				qstarts.clear();
				if (is_bed) { // bed specific
					std::vector<int> relative_blockstarts; // stores the block starts relative to the transcript start
					for (auto block : blockstarts) {
						relative_blockstarts.push_back(block - tstart);
					}

					bed_write_row(bed, prev_chrom, tstart, tend, qname, 
						1000, prev_strand, tstart, tend, "255,0,0", blockcount, 
						vector_to_str(blocksizes), vector_to_str(relative_blockstarts));
				} else { // psl specific
					int pos = 0;
					qstarts = {pos};

					for (int b = 0; b < blocksizes.size() - 1; b++) {
						pos += blocksizes[b];
						qstarts.push_back(pos);
					}

					// to get the query sequence size value for the pcl, we look in the dictionary
					int tsize = 0;
					if (chrom_to_size.size()) {
						if (chrom_to_size.count(prev_chrom)) {
							tsize = chrom_to_size[prev_chrom];
						} else {
							missing_chroms.insert(prev_chrom);
						}
					}

					psl_write_row(bed, 0, 0, 0, 0, 0, 0, 0, 0, prev_strand, qname, qsize, 0, qsize, prev_chrom, tsize, tstart, tend, (int)blockcount, vector_to_str(blocksizes), vector_to_str(qstarts), vector_to_str(blockstarts));
				
					}
			}

			// reset everything for the next transcript
			blockstarts.clear();
			blocksizes.clear();
			prev_transcript = this_transcript;
			// just get the gene_id from the attributes entry and assign it to prev_gene
			prev_gene = values[8].substr(values[8].find("gene_id") + 9, values[8].length());
			prev_gene = parseUntilChar(prev_gene, '"').first;
			prev_chrom = chrom;
			prev_strand = strand;
		}

		blockstarts.push_back(start);
		blocksizes.push_back(end - start);
  	}
	// last entry...
	blockcount = blockstarts.size();
	if (blockcount > 1 && blockstarts[0] > blockstarts[1]) { // need to reverse exons
		std::reverse(blocksizes.begin(), blocksizes.end());
		std::reverse(blockstarts.begin(), blockstarts.end());
	}

	// query (eg transcript)
	int qsize = std::accumulate(blocksizes.begin(), blocksizes.end(), 0);
	// target (eg chrom)
	tstart = blockstarts.front();
	tend = blockstarts.back() + blocksizes.back();
	std::string qname = this_transcript;


	if (is_bed) { // bed specific
		std::vector<int> relative_blockstarts; // stores the block starts relative to the transcript start
		for (auto block : blockstarts) {
			relative_blockstarts.push_back(block - tstart);
		}

		bed_write_row(bed, chrom, tstart, tend, qname, 
					1000, strand, start, end, "255,0,0", blockcount, 
					vector_to_str(blocksizes), vector_to_str(relative_blockstarts));

	}  else { // psl specific
		// to get the query sequence size value for the pcl, we look in the dictionary
		int tsize = 0;
		if (chrom_to_size.size()) {
			if (chrom_to_size.count(prev_chrom)) {
				tsize = chrom_to_size[prev_chrom];
			} else {
				missing_chroms.insert(prev_chrom);
			}
		}

		psl_write_row(bed, 0, 0, 0, 0, 0, 0, 0, 0, strand, qname, qsize, 0, qsize, chrom, tsize, tstart, tend, blockcount, vector_to_str(blocksizes), vector_to_str(qstarts), vector_to_str(blockstarts));
	}

	if (missing_chroms.size()) {
		Rcpp::Rcout << "chromosomes found in gtf but not in chrom_sizes file: " << set_to_str(missing_chroms) << "\n";
	}
	Rcpp::Rcout << "finished gtf_to_bed_cpp\n";

	gtf.close();
	bed.close();
	if (chrom_sizes.is_open()) chrom_sizes.close();
}