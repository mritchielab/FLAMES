#include "group_bam2isoform.h"

#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <sstream>
#include <algorithm>
#include <thread>

#include <Rcpp.h>

#include "../classes/GFFData.h"
#include "../classes/GeneAnnotationParser.h"
#include "../classes/GeneBlocks.h"
#include "../classes/Isoforms.h" 
#include "../classes/BamRecord.h"
#include "../utility/utility.h"
#include "../utility/bam.h"
#include "../utility/parsing.h"
#include "../utility/cigars.h"
// #include "../tests/test_utilities.h"

struct DataStruct2 {
	bam_header_t *header;
	Isoforms *isoform;
	int i;
};

#define BAM_CMATCH      0   // CIGAR character for matching
#define BAM_CDEL        2   
#define BAM_CREF_SKIP   3
#define BAM_CEQUAL      7
#define BAM_CDIFF       8

static std::unordered_map<std::string, std::string> 
get_fa_dict(std::string filename)
{
    // Modify in place
    auto
    string_toupper = [] (std::string &input)
    { std::transform(input.begin(), input.end(), input.begin(), [](unsigned char c) { return std::toupper(c); }); };

    std::string ch = "";
    std::stringstream sequence;
    std::unordered_map<std::string, std::string>  output;

    std::ifstream infile(filename);

    std::string line;
    while (std::getline(infile, line)) {
        if (line[0] == '>') {
            if (ch != "") {
                output[ch] = sequence.str();
            }
            ch = strip(std::string(line.begin()+1, line.end()));
            string_toupper(ch);
            sequence.str(""); // reset the stringstream to empty
        } else {
            string_toupper(line);
            sequence << strip(line);
        }
    }
    infile.close();

    return output;
}

static std::vector<StartEndPair>
get_blocks(const BAMRecord &record) {
    std::vector<StartEndPair> blocks;
    int pos = record.reference_start;

    for (const auto & [op, len] : record.cigar) {
        if ((op == BAM_CMATCH) ||
            (op == BAM_CEQUAL) ||
            (op == BAM_CDIFF)) {
            // add it to the blocks and move on
            blocks.push_back({pos, pos + len});
            pos += len;
        } else if (
            (op == BAM_CDEL) ||
            (op == BAM_CREF_SKIP)) {
            // just skip over this position
            pos += len;
        }
    }

    return blocks;
}

void 
merge_files_and_delete(std::ofstream &outfile, const std::vector<std::string> &infiles) {
    if (!outfile.is_open()) return;

    for (const auto &file : infiles) {
        std::ifstream infile(file);

        foreachLineinFile(file, [&outfile](const std::string &line) { 
            outfile << line << "\n";
        });

        std::remove(file.c_str());
    }
}
void merge_files_and_delete(const std::string &filename, const std::vector<std::string> &infiles) {
    std::ofstream outfile(filename);
    merge_files_and_delete(outfile, infiles);
    outfile.close();
}

static int
bam2isoform_fetch_function(const bam1_t *b, void *data)
{
    // retrieve the data from htslib and void * object
	DataStruct2 *data_struct = (DataStruct2 *)data;
	BAMRecord rec = read_record(b, data_struct->header);

    // smooth and update the CIGAR
	std::vector<CigarPair> cigar = smooth_cigar(rec.cigar, 20);
	rec.cigar = cigar;
	rec.cigar_string = generate_cigar(cigar);

	std::vector<StartEndPair> tmp_blocks = get_blocks(rec);
	Junctions junctions = blocks_to_junctions(tmp_blocks);
    
	std::vector<std::string> tmp_blocks_str = ranges::map<StartEndPair, std::string>(tmp_blocks, [](StartEndPair sep) -> std::string { return sep.getString(); });
    
    data_struct->isoform->add_isoform(junctions, rec.flag.read_reverse_strand);
	data_struct->i = data_struct->i + 1;

	return 0;
}

void group_bam2isoform(
    const std::string &bam_in,
    const std::string &out_gff3,
    const std::string &out_stat,
    const std::unordered_map<std::string, std::vector<StartEndPair>>  	&gene_dict, 
    const std::unordered_map<std::string, Junctions>           	 		&transcript_to_junctions,
    const std::unordered_map<std::string, Pos>                        	&transcript_dict,
    const std::unordered_map<std::string, std::vector<GeneBlocks>>      &chr_to_blocks,
    const std::string &fa_file,
    const Rcpp::List &isoform_parameters,
    const std::string &raw_gff3)
{
    // // random seed stuff here
    // if (false) {
    // // if (!file_exists(bam_in + ".bai")) { // WE CAN'T CHECK THIS AS THIS IS AN ABORT STATEMENT APPARENTLY
	// 	Rcpp::stop("Can not find corresponding .bai file %s. Cancelling group_bam2isoform.\n", bam_in);
    //     // Rcpp::Rcout << "Can not find corresponding .bai file " << bam_in << ". Cancelling group_bam2isoform.\n";
    //     return;
    // }
    
    // import all the values of fa_f
    const std::unordered_map<std::string, std::string> fa_dict = get_fa_dict(fa_file);

    std::string iso_annotated_output_prefix = getFilenameBeforeExt(out_gff3, '.');
    std::string tss_tes_output_prefix = getFilenameBeforeExt(out_stat, '.');
    std::string raw_gff3_prefix = raw_gff3 != "" ? getFilenameBeforeExt(raw_gff3, '.') : "";
    std::vector<std::thread> pool;
    std::vector<std::string> isoform_annotated_files;
    std::vector<std::string> tss_tes_files;
    std::vector<std::string> raw_gff3_files;

    // capture variables form Rcpp::List isoform_parameters
    int max_ts_dist=isoform_parameters["max_ts_dist"],
        max_dist=isoform_parameters["max_dist"],
        strand_specific=isoform_parameters["strand_specific"],
        max_splice_match_dist=isoform_parameters["max_splice_match_dist"],
        remove_incomp_reads=isoform_parameters["remove_incomp_reads"],
        min_fl_exon_len=isoform_parameters["min_fl_exon_len"],
        min_sup_cnt=isoform_parameters["min_sup_cnt"], 
        max_site_per_splice=isoform_parameters["max_site_per_splice"];
        
    double min_sup_pct=isoform_parameters["min_sup_pct"],
        min_cnt_pct=isoform_parameters["min_cnt_pct"];

    // run the main group_bam2isoform on each chromosome,
    // generating a new worker thread to handle the isoform identification
    // and merging the resultant out files into one output gff3
    for (const auto &[chromosome, blocks] : chr_to_blocks) {
        std::string isoform_outfile = iso_annotated_output_prefix + chromosome + ".gff3";
        isoform_annotated_files.push_back(isoform_outfile);
        std::string tss_tes_outfile = tss_tes_output_prefix + chromosome + ".bedgraph";
        tss_tes_files.push_back(tss_tes_outfile);
        std::string raw_gff3_outfile = raw_gff3 != "" ? raw_gff3_prefix + chromosome + ".gff3" : "";
        raw_gff3_files.push_back(raw_gff3_outfile);

        auto thread_function = [
            chromosome=chromosome, &blocks=blocks, 
            bam_in,
            &transcript_to_junctions, 
            &transcript_dict, &gene_dict, &fa_dict,
            raw_gff3_outfile, tss_tes_outfile, isoform_outfile,
            max_ts_dist,
            max_dist,
            strand_specific,
            max_splice_match_dist,
            remove_incomp_reads,
            min_fl_exon_len,
            min_sup_cnt, 
            min_sup_pct,
            max_site_per_splice,
            min_cnt_pct
        ]() {
            bamFile bam = bam_open(bam_in.c_str(), "r"); // bam.h
            bam_index_t *bam_index = bam_index_load(bam_in.c_str());
            bam_header_t *header = bam_header_read(bam); // bam.h
            int tid = bam_get_tid(header, chromosome.c_str());
            
            for (const auto &block : blocks) {
                Isoforms tmp_isoform(
                    chromosome, max_ts_dist,
                    max_dist,
                    strand_specific,
                    max_splice_match_dist,
                    remove_incomp_reads,
                    min_fl_exon_len,
                    min_sup_cnt, min_sup_pct,
                    max_site_per_splice);
                DataStruct2 data = {header, &tmp_isoform, 0};

                bam_fetch(bam, bam_index, tid, block.start, block.end, &data, &bam2isoform_fetch_function);
                
                auto TSS_TES_site = get_TSS_TES_site(transcript_to_junctions, block.transcript_list);

                // then process the complete isoform
                if (tmp_isoform.size() > 0) {
                    std::ofstream tss_tes_stream(tss_tes_outfile);

                    tmp_isoform.update_all_splice();
                    tmp_isoform.filter_TSS_TES(tss_tes_stream, TSS_TES_site, (float)0.1);
                    
                    tmp_isoform.match_known_annotation(
                        transcript_to_junctions,
                        transcript_dict,
                        gene_dict,
                        block,
                        fa_dict
                    );
                    
                    if (raw_gff3_outfile != "") {
                        // todo - i haven't written the Isoforms function to do this yet
                        // splice_raw.write(tmp_isoform()); 
                    }

                    tss_tes_stream.close();

                    std::ofstream iso_annotated(isoform_outfile);
                    iso_annotated << tmp_isoform.isoform_to_gff3(min_cnt_pct);
                    iso_annotated.close();
                }
            }

            bam_close(bam);
        }; // end thread_function

        pool.push_back(std::thread(thread_function));
    }

    for (auto &thread : pool) {
        thread.join();
    }

    // After joining all threads 
    // merge all the individual output files into a single output, and delete all temp files
    std::ofstream iso_annotated(out_gff3);
    iso_annotated << "##gff-version 3\n";
    merge_files_and_delete(iso_annotated, isoform_annotated_files);
    iso_annotated.close();

    merge_files_and_delete(out_stat, tss_tes_files);
    // add to splice_raw if we are planning on outputting raw_gff3
    std::ofstream splice_raw;
    if (raw_gff3 != "") {
        splice_raw.open(raw_gff3);
        splice_raw << "##gff-version 3\n";
        merge_files_and_delete(splice_raw, raw_gff3_files);
        splice_raw.close();
    }
}