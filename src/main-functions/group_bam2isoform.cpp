#include "group_bam2isoform.h"

#include <string>
#include <map>
#include <vector>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <Rcpp.h>

#include "config.h"
#include "GeneBlocks.h"
#include "junctions.h"
#include "misc.h"
#include "Isoforms.h"
#include "bam.h"
#include "StartEndPair.hpp"
#include "cigars.h"
#include "BamRecord.hpp"

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

// [[Rcpp::export]]
void
bam_read (std::string bam_in, std::string chr, int s, int e)
{
    // read a bamfile
    bamFile bam = bam_open(bam_in.c_str(), "r"); // bam.h
    bam_index_t *bam_index = bam_index_load(bam_in.c_str());
    bam_header_t *header = bam_header_read(bam); // bam.h
    int tid = bam_get_tid(header, chr.c_str());

    

    Rcpp::Rcout << "bamfile is done\n";
    Rcpp::Rcout << "tid is " << tid << "\n";

    std::vector<BAMRecord>
    records;
    DataStruct data = {header, &records};
    auto it_region = bam_fetch(bam, bam_index, tid, s, e, &data, fetch_function);

    bam_close(bam);
}

std::vector<StartEndPair>
get_blocks(BAMRecord record)
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
            pos += 1;
        } else if (
            (op == BAM_CDEL) ||
            (op == BAM_CREF_SKIP)) {
            // just skip over this position
            pos += 1;
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

            auto it = bam_fetch(bam, bam_index, tid, block.start, block.end, &data, &fetch_function);
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
    std::unordered_map<std::string, std::vector<GeneBlocks>>    * chr_to_blocks, 
    std::unordered_map<std::string, std::vector<StartEndPair>>  * gene_dict, 
    std::unordered_map<std::string, Junctions>                  * transcript_to_junctions,
    std::unordered_map<std::string, Pos>                        * transcript_dict,
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

    // add to splice_raw if we are planning on outputting raw_gff3
    std::ofstream splice_raw;
    if (raw_gff3 != "") {
        splice_raw.open(raw_gff3);
        splice_raw << "##gff-version 3\n";
    }

    std::ofstream 
    tss_tes_stat (out_stat);

    std::unordered_map<IsoformKey, Isoforms*>
    isoform_dict = {};

    // import all the values of fa_f
    std::unordered_map<std::string, std::string>
    fa_dict;
    for (const auto & c : get_fa(fa_f)) {
        fa_dict[c.first] = c.second;
    }

    std::vector<BAMRecord>
    records = {};
    DataStruct data = {header, &records};

    Rcpp::Rcout << "\tlooking at chr_to_blocks of len " << chr_to_blocks->size() << "\n";
    for (const auto & [chr, blocks] : (*chr_to_blocks)) {
        int tid = bam_get_tid(header, chr.c_str());
        Rcpp::Rcout << "\t\tstarting on " << chr << " which has " << blocks.size() << " blocks\n";

        int ith = 0;
        for (const auto & block : blocks) {
            Rcpp::Rcout << "looking at " << ith << " block: (" << block.start << ", " << block.end << ")\n";
            ith++;

            records = {};
            // extract this from the bam file

            auto it = bam_fetch(bam, bam_index, tid, block.start, block.end, &data, &fetch_function);
            auto TSS_TES_site = get_TSS_TES_site(transcript_to_junctions, &(block.transcript_list));
            auto tmp_isoform = new Isoforms (chr, isoform_parameters);
            
            // add all the records in the bamfile to the Isoform object
            int recnum = 0;
            for (const auto & rec : records) {
                Rcpp::Rcout << "\tlooking at rec " << recnum << "\n";
                recnum++;
				
                std::vector<CigarPair>
                cigar_smooth = smooth_cigar(rec.cigar, 20);
                std::string
                cigar_string = generate_cigar(rec.cigar);
                std::vector<StartEndPair>
                tmp_blocks = get_blocks(rec);
                Junctions
                junctions = blocks_to_junctions(tmp_blocks);
                
                Rcpp::Rcout << "\t\tadding junctions to isoform:\n{'right':" << junctions.right; 
                Rcpp::Rcout << ", 'junctions': (";
                for (const auto & j : junctions.junctions) {
                    Rcpp::Rcout << j << ", ";
                }
                Rcpp::Rcout << "), 'left': " << junctions.left << "}\n";
                tmp_isoform->add_isoform(junctions, rec.flag.read_reverse_strand);
            }

            // then process the isoform
            if (tmp_isoform->len() > 0) {
                Rcpp::Rcout << "\t\tfound an isoform with len>0\n";
                tmp_isoform->update_all_splice();
                Rcpp::Rcout << "finished update_all_splice\n";
                tmp_isoform->filter_TSS_TES(&tss_tes_stat, TSS_TES_site, (float)0.1);
                Rcpp::Rcout << "finished filter_TSS_TES\n";
                tmp_isoform->match_known_annotation(
                    (*transcript_to_junctions),
                    (*transcript_dict),
                    (*gene_dict),
                    block,
                    fa_dict
                );
                Rcpp::Rcout << "made it to isoform_dict adding\n";
                // isoform_dict[{chr, block.start, block.end}] = Isoforms(chr, isoform_parameters);
                isoform_dict[{chr, block.start, block.end}] = tmp_isoform;
                
                if (raw_gff3 != "") {
                    // todo - i haven't written the Isoforms function to do this yet
                    // splice_raw.write(tmp_isoform()); 
                }
                Rcpp::Rcout << "made it to isoform_annotated adding\n";
                iso_annotated << tmp_isoform->isoform_to_gff3(isoform_parameters.MIN_CNT_PCT);
            }
            Rcpp::Rcout << "made it to deletion\n";
            delete tmp_isoform;
        }
    }

    // finally, close all the files
    bam_close(bam);
    tss_tes_stat.close();
    iso_annotated.close();
    if (raw_gff3 != "") {
        splice_raw.close();
    }
}