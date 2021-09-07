void
group_bam2isoform (
    std::string bam_in, 
    std::string out_gff3, 
    std::string out_stat, 
    std::map<std::string, 
    std::vector<GeneBlocks>> chr_to_blocks, 
    std::map<std::string, std::vector<std::pair<int, int>>> gene_dict, 
    std::map<std::string, Junctions> transcript_to_junctions,
    std::map<std::string, Pos> transcript_dict,
    std::string fa_f,
    std::map<std::string, int> config,
    std::string raw_gff3 = ""
)
{
    if (config.count("random_seed"))
    {
        srand(config["random_seed"]);
    }
    else
    {
        srand(666666);
    }

    // i need to learn to read a bamfile


    // set up all the output files
    std::ofstream iso_annotated;
    iso_annotated.open(out_gff3);
    iso_annotated << "##gff-version 3\n";

    if (raw_gff3 != "")
    {
        std::ofstream splice_raw;
        splice_raw.open(raw_gff3);
        splice_raw << "##gff-version 3\n";
    }

    std::ofstream tss_tes_stat;
    tss_tes_stat.open(out_stat);

    std::map<std::vector<int>, Iso>
    isoform_dict;

    // import all the values of fa_f
    std::map<std::string, std::string>
    fa_dict;
    for (const auto & c : get_fa(fa_f))
    {
        fa_dict[c.first] = c.second;
    }

    for (const auto & [chr, blocks] : chr_to_blocks)
    {
        for (const auto & block : blocks)
        {
            // extract this from the bam file
            auto it_region = 0;

            auto TSS_TES_site = get_TSS_TES_site(transcript_to_junctions, block.transcript_list);
            auto tmp_isoform = Isoforms(chr, config);

        }
    }
}