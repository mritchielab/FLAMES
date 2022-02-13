#include "group_bam2isoform.h"

#include <string>
#include <map>
#include <vector>
#include <fstream>
#include <Rcpp.h>

#include "../classes/Config.h"
#include "../classes/GeneBlocks.h"
#include "../classes/Isoforms.h"
#include "../classes/StartEndPair.h"
#include "../classes/BamRecord.h"
#include "../classes/IsoKey.h"
#include "../classes/DataStruct.h"
#include "../utility/junctions.h"
#include "../utility/misc.h"
#include "../utility/bam.h"
#include "../utility/cigars.h"



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

// void
// bam_read (std::string bam_in, std::string chr, int s, int e)
// {
//     // read a bamfile
//     bamFile bam = bam_open(bam_in.c_str(), "r"); // bam.h
//     bam_index_t *bam_index = bam_index_load(bam_in.c_str());
//     bam_header_t *header = bam_header_read(bam); // bam.h
//     int tid = bam_get_tid(header, chr.c_str());
 	
    
// 	// the goal is to read in a BAM file and acquire a vector of BAMRecord
// 	// of every BAM entry in the BAM file
// 	std::vector<BAMRecord> records;
// 	DataStruct data = {header, &records};
//     // auto it_region = bam_fetch(bam, bam_index, tid, s, e, &data, fetch_function);
// 	bam_fetch(bam, bam_index, tid, s , e, &data, fetch_function);
//     bam_close(bam);
// }

std::vector<StartEndPair>
get_blocks(const BAMRecord &record)
{
    std::vector<StartEndPair>
    blocks = {};
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
minimal_group_bam2isoform
(

    std::string bam_in, 
    std::string out_gff3, 
    std::string out_stat, 
    std::unordered_map<std::string, std::vector<GeneBlocks>>    * chr_to_blocks, 
    std::unordered_map<std::string, std::vector<StartEndPair>>  * gene_dict, 
    std::unordered_map<std::string, Junctions>                  * transcript_to_junctions,
    std::unordered_map<std::string, Pos>                        * transcript_dict,
    std::string fa_f,
    IsoformParameters isoform_parameters,
    std::string raw_gff3
)
{
    // read a bamfile
    bamFile bam = bam_open(bam_in.c_str(), "r"); // bam.h
    bam_index_t *bam_index = bam_index_load(bam_in.c_str());
    bam_header_t *header = bam_header_read(bam); // bam.h


    std::vector<BAMRecord>
    records = {};
    Rcpp::Rcout << "\tmade records\n";
    DataStruct data = {header, &records};


    Rcpp::Rcout << "made it to the big for\n";
    for (const auto & [chr, blocks] : (*chr_to_blocks)) {
        Rcpp::Rcout << "started loop with " << chr << "\n";
        int tid = bam_get_tid(header, chr.c_str());

        for (const auto & block : blocks) {
            Rcpp::Rcout << "started inner loop with " << block.start << ", " << block.end << "\n";
            // if (block.start == 8208472) {
            //     Rcpp::Rcout << "!! SKIPPING the danger zone !!\n";
            //     continue;
            // }
            Rcpp::Rcout << "\tcleared records";
            // extract this from the bam file

            Rcpp::Rcout << "\tabout to bamfetch\n";

            // auto it = bam_fetch(bam, bam_index, tid, block.start, block.end, &data, &fetch_function);
            bam_fetch(bam, bam_index, tid, block.start, block.end, &data, &fetch_function);
			Rcpp::Rcout << "bam_fetch done\n";
        }
    }


    bam_close(bam);
}

void
group_bam2isoform (
    std::string bam_in, 
    std::string out_gff3, 
    std::string out_stat, 
    const std::unordered_map<std::string, std::vector<GeneBlocks>> 		& chr_to_blocks, 
    const std::unordered_map<std::string, std::vector<StartEndPair>>  	& gene_dict, 
    const std::unordered_map<std::string, Junctions>           	 		& transcript_to_junctions,
    const std::unordered_map<std::string, Pos>                        	& transcript_dict,
    std::string fa_f,
    IsoformParameters isoform_parameters,
    std::string raw_gff3
)
{
    // if (config.count("random_seed")) {
    //     srand(config["random_seed"]);
    // } else {
    //     srand(666666);
    // }

    // read a bamfile
    bamFile bam = bam_open(bam_in.c_str(), "r"); // bam.h
    bam_index_t *bam_index = bam_index_load(bam_in.c_str());
    bam_header_t *header = bam_header_read(bam); // bam.h

    // set up all the output files
    std::ofstream iso_annotated;
    iso_annotated.open(out_gff3);
    iso_annotated << "##gff-version 3\n";

    std::ofstream also_out ("also_out_cpp.gff");
    also_out << "##gff-version 3\n";
    // add to splice_raw if we are planning on outputting raw_gff3
    std::ofstream splice_raw;
    if (raw_gff3 != "") {
        splice_raw.open(raw_gff3);
        splice_raw << "##gff-version 3\n";
    }

    std::ofstream 
    tss_tes_stat (out_stat);

    // std::unordered_map<IsoformKey, Isoforms*>
    // isoform_dict = {};

    // import all the values of fa_f
    std::unordered_map<std::string, std::string>
    fa_dict;
    for (const auto & c : get_fa(fa_f)) {
        fa_dict[c.first] = c.second;
    }

    std::vector<BAMRecord>
    records = {};
    DataStruct data = {header, &records};

    for (const auto & [chr, blocks] : chr_to_blocks) {
        int tid = bam_get_tid(header, chr.c_str());


        int ith = 0;
        for (const auto & block : blocks) {
            ith++;
            records = {};
            // extract this from the bam file

            // auto it = bam_fetch(bam, bam_index, tid, block.start, block.end, &data, &fetch_function);
            bam_fetch(bam, bam_index, tid, block.start, block.end, &data, &fetch_function);
			auto TSS_TES_site = get_TSS_TES_site(transcript_to_junctions, block.transcript_list);
            Isoforms tmp_isoform(chr, isoform_parameters);
            
            // add all the records in the bamfile to the Isoform object
            int recnum = 0;

            for (BAMRecord & rec : records) {
                recnum++;
				
                auto cigar = smooth_cigar(rec.cigar, 20);
                rec.cigar = cigar;
                std::string
                cigar_string = generate_cigar(cigar);
                std::vector<StartEndPair>
                tmp_blocks = get_blocks(rec);
                Junctions
                junctions = blocks_to_junctions(tmp_blocks);
                
                tmp_isoform.add_isoform(junctions, rec.flag.read_reverse_strand);
            }

            // then process the isoform
            if (tmp_isoform.size() > 0) {
                tmp_isoform.update_all_splice();
                tmp_isoform.filter_TSS_TES(tss_tes_stat, TSS_TES_site, (float)0.1);
                tmp_isoform.match_known_annotation(
                    transcript_to_junctions,
                    transcript_dict,
                    gene_dict,
                    block,
                    fa_dict
                );
                // isoform_dict[{chr, block.start, block.end}] = Isoforms(chr, isoform_parameters);
                // isoform_dict[{chr, block.start, block.end}] = tmp_isoform;
                
                if (raw_gff3 != "") {
                    // todo - i haven't written the Isoforms function to do this yet
                    // splice_raw.write(tmp_isoform()); 
                }
                iso_annotated << tmp_isoform.isoform_to_gff3(isoform_parameters.MIN_CNT_PCT);
            }
        }
    }

    // finally, close all the files
    bam_close(bam);
    tss_tes_stat.close();
    iso_annotated.close();
    also_out.close();
    if (raw_gff3 != "") {
        splice_raw.close();
    }
}