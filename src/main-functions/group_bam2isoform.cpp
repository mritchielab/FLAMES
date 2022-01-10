#include "group_bam2isoform.h"

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

    std::cout << "bamfile is done\n";
    std::cout << "tid is " << tid << "\n";

    std::vector<BAMRecord>
    records;
    DataStruct data = {header, &records};
    auto it_region = bam_fetch(bam, bam_index, tid, s, e, &data, fetch_function);
    std::cout << "found " << records.size() << " records\n";
    bam_close(bam);
}

std::vector<StartEndPair>
get_blocks(int reference_start, std::vector<CigarPair> cigar)
{
    std::vector<StartEndPair>
    blocks = {};

    int pos = reference_start;

    for (const auto & [op, len] : cigar) {
        if ((op == BAM_CMATCH) ||
            (op == BAM_CEQUAL) ||
            (op == BAM_CDIFF )) {
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
    std::cout << "\tmade records\n";
    DataStruct data = {header, &records};


    std::cout << "made it to the big for\n";
    for (const auto & [chr, blocks] : (*chr_to_blocks)) {
        std::cout << "started loop with " << chr << "\n";
        int tid = bam_get_tid(header, chr.c_str());

        for (const auto & block : blocks) {
            std::cout << "started inner loop with " << block.start << ", " << block.end << " (" << "\n";
            // if (block.start == 8208472) {
            //     std::cout << "!! SKIPPING the danger zone !!\n";
            //     continue;
            // }
            std::cout << "\tcleared records";
            // extract this from the bam file

            std::cout << "\tabout to bamfetch\n";

            auto it = bam_fetch(bam, bam_index, tid, block.start, block.end, &data, &fetch_function);
            std::cout << "bam_fetch done\n";
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
    std::ofstream
    also_out_gff3 ("also_out.gff");
    also_out_gff3 << "##gff-version 3\n";

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
    std::cout << "fa_dict.size():" << fa_dict.size() << "\n";

    std::vector<BAMRecord>
    records = {};
    DataStruct data = {header, &records};

    std::cout << "\tlooking at chr_to_blocks of len " << chr_to_blocks->size() << "\n";
    for (const auto & [chr, blocks] : (*chr_to_blocks)) {
        int tid = bam_get_tid(header, chr.c_str());
        std::cout << "\t\tstarting on " << chr << " which has " << blocks.size() << " blocks\n";

        int ith = 0;
        for (const auto & block : blocks) {
            std::cout << "looking at " << ith << " block: (" << block.start << ", " << block.end << ")\n";
            std::cout << "\tblock.transcript_list (size " << block.transcript_list.size() << "):";
            for (const auto & tr : block.transcript_list) {
                std::cout << tr << ",";
            }
            std::cout << "\n";
            ith++;

            records = {};
            // extract this from the bam file
            auto it = bam_fetch(bam, bam_index, tid, block.start, block.end, &data, &fetch_function);
            std::cout << "found " << records.size() << " records\n";
            auto TSS_TES_site = get_TSS_TES_site(transcript_to_junctions, &(block.transcript_list));
            std::cout << "TSS_TES_site=";
            junctions_print(TSS_TES_site);

            auto tmp_isoform = new Isoforms (chr, isoform_parameters);
            
            // add all the records in the bamfile to the Isoform object
            int recnum = 0;
            for (const auto & rec : records) {
                std::cout << "\tlooking at rec " << recnum << "\n";
                recnum++;

                std::cout << "\t\t\tprior to smoothing, cigar is " << rec.cigar.size() << "\n";
                auto
                cigar_smooth = smooth_cigar(rec.cigar, 20);

                std::cout << "\t\t\tlen(cigar) is " << cigar_smooth.size() << ", cigar is [";
                for (const auto & c : cigar_smooth) {
                    std::cout << "[" << c.op << "," << c.len << "] ";
                }
                std::cout << "]\n";

                auto
                cigar_string = generate_cigar(cigar_smooth);
                std::cout << "\t\t\tcigarstring is " << cigar_string << "\n";
                auto
                tmp_blocks = get_blocks(rec.reference_start, cigar_smooth);
                std::cout << "\t\t\ttmp_blocks size is " << tmp_blocks.size() << "\n";
                auto
                junctions = blocks_to_junctions(tmp_blocks);
                
                std::cout << "\t\tadding junctions to isoform:\n";
                junctions_print(junctions);
                
                tmp_isoform->add_isoform(junctions, rec.flag.read_reverse_strand);
            }

            // then process the isoform
            std::cout << "\t\tup to the if\n";
            if (tmp_isoform->len() > 0) {
                std::cout << "\t\tfound an isoform with len>0\n";
                tmp_isoform->update_all_splice();
                std::cout << "finished update_all_splice\n";
                tmp_isoform->filter_TSS_TES(&tss_tes_stat, TSS_TES_site, (float)0.1);
                std::cout << "finished filter_TSS_TES\n";
                tmp_isoform->match_known_annotation(
                    (*transcript_to_junctions),
                    (*transcript_dict),
                    (*gene_dict),
                    block,
                    fa_dict
                );
                std::cout << "\t\t\tstop! isoform check.\n";
                std::cout << "\t\t\t\tjunction_dict.size(): " << tmp_isoform->junction_dict.size() << "\n";
                std::cout << "\t\t\t\tjunction_list.size(): " << tmp_isoform->junction_list.size() << "\n";
                std::cout << "\t\t\t\tlr_pair.size(): " << tmp_isoform->lr_pair.size() << "\n";
                std::cout << "\t\t\t\tleft.size(): " << tmp_isoform->left.size() << "\n";
                std::cout << "\t\t\t\tright.size(): " << tmp_isoform->right.size() << "\n";
                std::cout << "\t\t\t\tsingle_block_dict.size(): " << tmp_isoform->single_block_dict.size() << "\n";
                std::cout << "\t\t\t\tsingle_blocks.size(): " << tmp_isoform->single_blocks.size() << "\n";
                std::cout << "\t\t\t\tstrand_counts.size(): " << tmp_isoform->strand_counts.size() << "\n";
                std::cout << "\t\t\t\tknown_isoforms.size(): " << tmp_isoform->known_isoforms.size() << "\n";
                std::cout << "\t\t\t\tnew_isoforms.size(): " << tmp_isoform->new_isoforms.size() << "\n";
                std::cout << "\t\t\t\traw_isoforms.size(): " << tmp_isoform->raw_isoforms.size() << "\n";
                std::cout << "\t\t\t\tge_dict.size(): " << tmp_isoform->ge_dict.size() << "\n";
                std::cout << "made it to isoform_dict adding\n";
                // isoform_dict[{chr, block.start, block.end}] = Isoforms(chr, isoform_parameters);
                isoform_dict[{chr, block.start, block.end}] = tmp_isoform;
                
                if (raw_gff3 != "") {
                    // todo - i haven't written the Isoforms function to do this yet
                    // splice_raw.write(tmp_isoform()); 
                }
                std::cout << "made it to isoform_annotated adding\n";
                iso_annotated << tmp_isoform->isoform_to_gff3(isoform_parameters.MIN_CNT_PCT);
                also_out_gff3 << tmp_isoform->isoform_to_gff3(isoform_parameters.MIN_CNT_PCT);
            }
            std::cout << "made it to deletion\n";
            delete tmp_isoform;
        }
    }
    std::cout << "finished group_bam2isoform, isoform_dict.size() is " << isoform_dict.size() << "\n";
    // finally, close all the files
    bam_close(bam);
    tss_tes_stat.close();
    iso_annotated.close();
    also_out_gff3.close();
    if (raw_gff3 != "") {
        splice_raw.close();
    }
}