#include "group_bam2isoform.h"

static int
fetch_function(const bam1_t *b, void *data)
{
    int * found = (int*)data;
    // print out some output
    std::cout << b->data << ", " << b->core.tid << ", " << b->core.pos << ", " << b->core.l_qseq << "\n";
    // and change i to 6 to signify that we've found something
    (*found)++;
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

    std::cout << "bamfile is done\n";
    
    int found = 0;

    std::cout << "found is " << found << "\n";
    auto it_region = bam_fetch(bam, bam_index, 0, s, e, &found, fetch_function);
    std::cout << "now found is " << found << " and it region is " << it_region << "\n";

    bam_close(bam);
}

void
group_bam2isoform (
    std::string bam_in, 
    std::string out_gff3, 
    std::string out_stat, 
    std::map<std::string, std::vector<GeneBlocks>>          * chr_to_blocks, 
    std::map<std::string, std::vector<std::pair<int, int>>> * gene_dict, 
    std::map<std::string, Junctions>                        * transcript_to_junctions,
    std::map<std::string, Pos>                              * transcript_dict,
    std::string fa_f,
    IsoformParameters isoform_parameters,
    std::string raw_gff3 = ""
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
            auto it_region = bam_fetch(bam, bam_index, 0, block.start, block.end, 0, fetch_function);

            auto TSS_TES_site = get_TSS_TES_site((*transcript_to_junctions), block.transcript_list);
            auto tmp_isoform = Isoforms(chr, isoform_parameters);
        }
    }
}