#include "group_bam2isoform.h"

/*  take a bam entry,
    populate a vector with all of its CIGAR operations
    this is to mimic the output of bamnostic's default cigar from bamfile.fetch
*/
std::vector<CigarPair>
generate_cigar_pairs(const bam1_t *b)
{
    std::vector<CigarPair>
    cigar_pairs;

    // iterate over the cigar
    const auto cigar = bam_get_cigar(b);
    for (int k = 0; k < b->core.n_cigar; k++) {
        cigar_pairs.push_back((CigarPair){
            bam_cigar_op(cigar[k]),
            bam_cigar_oplen(cigar[k])
        });
    }
    return cigar_pairs;
}

static int
fetch_function(const bam1_t *b, void *data)
{
    Record
    rec;

    std::vector<Record> (*records) = (std::vector<Record>*)data;

    std::cout << b->data << ", " << b->core.tid << ", " << b->core.pos << ", " << b->core.n_cigar << "\n";

    rec.reference_start = b->core.pos;
    rec.reference_end = b->core.pos + b->core.l_qseq;
    rec.is_reverse = b->core.flag;
    rec.cigar = generate_cigar_pairs(b);

    records->push_back(rec);
	return 0;
}

// [[Rcpp::export]]
void
bam_read (std::string bam_in, int s, int e)
{
    // read a bamfile
    bamFile bam = bam_open(bam_in.c_str(), "r"); // bam.h
    bam_index_t *bam_index = bam_index_load(bam_in.c_str());
    bam_header_t *header = bam_header_read(bam); // bam.h
    int tid = bam_get_tid(header, "chr22");

    std::cout << "bamfile is done\n";
    std::cout << "tid is " << tid << "\n";

    std::vector<Record>
    records;
    auto it_region = bam_fetch(bam, bam_index, tid, s, e, &records, fetch_function);

    bam_close(bam);
}

std::vector<StartEndPair>
get_blocks(Record record)
{
    std::vector<StartEndPair>
    blocks = {};

    int pos = record.reference_start;

    for (const auto & [op, len] : record.cigar) {
        if ((op == BAM_CMATCH) ||
            (op == BAM_CEQUAL) ||
            (op == BAM_CDIFF)) {
            // add it to the blocks and move on
            blocks.push_back({pos, pos + len})
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
    std::string raw_gff3 =
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
    bam_close(bam);

    // set up all the output files
    std::ofstream iso_annotated;
    iso_annotated.open(out_gff3);
    iso_annotated << "##gff-version 3\n";

    // add to splice_raw if we are planning on outputting raw_gff3
    if (raw_gff3 != "") {
        std::ofstream splice_raw;
        splice_raw.open(raw_gff3);
        splice_raw << "##gff-version 3\n";
    }

    std::ofstream 
    tss_tes_stat (out_stat);

    std::map<std::vector<int>, Iso>
    isoform_dict;

    // import all the values of fa_f
    std::map<std::string, std::string>
    fa_dict;
    for (const auto & c : get_fa(fa_f)) {
        fa_dict[c.first] = c.second;
    }

    for (const auto & [chr, blocks] : (*chr_to_blocks)) {
        for (const auto & block : blocks) {
            // extract this from the bam file

            int tid = bam_get_tid(header, chr.c_str());
            
            std::vector<Record>
            records = {};

            auto it_region = bam_fetch(bam, bam_index, tid, block.start, block.end, &records, fetch_function);
            auto TSS_TES_site = get_TSS_TES_site((*transcript_to_junctions), block.transcript_list);
            auto tmp_isoform = Isoforms(chr, isoform_parameters);

            // add all the records in the bamfile to the Isoform object
            for (const auto & rec : records) {
                auto
                cigar_smooth = smooth_cigar(rec.cigar, 20);
                auto
                cigar_string = generate_cigar(rec.cigar);
                auto
                tmp_blocks = get_blocks(rec);
                auto
                junctions = blocks_to_junctions(tmp_blocks);
                tmp_isoform.add_isoform(junctions, rec.is_reverse);
            }

            // then process the isoform
            if (tmp_isoform.len() > 0) {
                tmp_isoform.update_all_splice();
                tmp_isoform.filter_TSS_TES(tss_tes_stat, TSS_TES_site, 0.1);
                tmp_isoform.match_known_annotation(
                    (*transcript_to_junctions),
                    (*transcript_dict),
                    (*gene_dict),
                    block,
                    fa_dict
                );
                // isoform_dict[]
            }
        }
    }
}