#include "find_isoform.h"

Rcpp::List
isoform_objects_to_R(IsoformObjects * isoform_objects)
{
    auto wrap_dict = [] (std::unordered_map<std::string, Pos> dict)
    {
        Rcpp::CharacterVector keys;
        Rcpp::List list;

        for (auto & [key, pos] : dict) {
            keys.push_back(key);
            list.push_back(pos_to_R(&pos));
        }

        list.attr("names") = keys;

        return list;
    };

    return Rcpp::List::create(
        Rcpp::Named("transcript_dict") = wrap_dict(isoform_objects->transcript_dict),
        Rcpp::Named("transcript_dict_iso") = wrap_dict(isoform_objects->transcript_dict_iso)
    );
}

IsoformObjects
isoform_objects_from_R(Rcpp::List list)
{
    IsoformObjects isoform_objects;

    auto unwrap_dict = [] (Rcpp::List list)
    {
        std::unordered_map<std::string, Pos>
        dict;

        Rcpp::CharacterVector keys = list.attr("names");

        for (int i = 0; i < list.size(); ++i) {
            dict[(char*)keys[i]] = pos_from_R(list[i]);
        }

        return dict;
    };

    isoform_objects.transcript_dict = unwrap_dict(list["transcript_dict"]);
    isoform_objects.transcript_dict_iso = unwrap_dict(list["transcript_dict_iso"]);

    return isoform_objects;
}

void
log_params
(
    std::unordered_map<std::string, std::vector<GeneBlocks>> *
    chr_to_blocks,

    std::unordered_map<std::string, std::vector<StartEndPair>> *
    gene_dict,

    std::unordered_map<std::string, Junctions> *
    transcript_to_junctions,

    std::unordered_map<std::string, Pos> *
    transcript_dict
)
{
    std::ofstream
    outfile ("outfile.txt");

    outfile << "chr_to_blocks:\n";
    for (const auto & [chr, blocks] : *chr_to_blocks) {
        outfile << "\t" << chr << ":\n";
        for (const auto & block : blocks) {
            outfile << "\t\t(" << block.start << "," << block.end << ")\n";
        }
    }

    outfile << "\ngene_dict:\n";
    for (const auto & [gene, pairs] : *gene_dict) {
        outfile << "\t" << gene << ":\n";
        std::cout << "doing gene " << gene << "\n";
        if (pairs.size() == 0) {
            std::cout << "no pairs\n";
        }
        for (const auto & pair : pairs) {
            outfile << "\t\t(" << pair.start << "," << pair.end << ")\n";
        }
    }

    outfile << "\ntranscript_to_junctions:\n";
    for (auto & [tr, junctions] : *transcript_to_junctions) {
        outfile << "\t" << tr << ":{'right': " << junctions.right.front() << ", 'junctions': (";
        if (junctions.junctions.size() > 0) {
            for (auto it = junctions.junctions.begin(); it != junctions.junctions.end() - 1; ++it) {
                outfile << *it << ", ";
            }
            outfile << junctions.junctions.back();
        }
        outfile << "), 'left': " << junctions.left.front() << "}\n";
    }

    outfile << "\ntranscript_dict:\n";
    for (const auto & [tr, pos] : *transcript_dict) {
        outfile << "\t" << tr << ":\n\t\t(" << pos.start << "," << pos.end << ")\n";
    }
}

// [[Rcpp::export]]
Rcpp::List
find_isoform_cpp
(
    std::string gff3, 
    std::string genome_bam, 
    std::string isoform_gff3, 
    std::string tss_test_stat, 
    std::string genomefa, 
    std::string transcript_fa, 
    int         downsample_ratio, 
    Rcpp::List  config_list, 
    std::string raw_splice_isoform
)
{
    std::cout << "#### Reading Gene Annotations\n";

    // first, extract from the gff file
    GFFData gene_anno = parse_gff_or_gtf(gff3);
    log_gff_data(gene_anno, "initial_gff_log_cpp.txt");

    auto chr_to_gene        = gene_anno.chr_to_gene;
    std::cout << "chr_to_gene is size " << chr_to_gene.size() << "\n";
    for (const auto & [chr, gene] : chr_to_gene) {
        std::cout << "\t" << chr << ":" << gene.size() << "\n";
    }

    std::cout << "chr_to_gene is currently " << chr_to_gene.size() << " long\n";
    auto transcript_dict    = gene_anno.transcript_dict;
    std::cout << "transcript_dict is currently " << transcript_dict.size() << " long\n";
    auto gene_to_transcript = gene_anno.gene_to_transcript;
    std::cout << "gene_to_transcript is currently " << gene_to_transcript.size() << " long\n";
    for (const auto & [gene, transcript] : gene_to_transcript) {
        std::cout << "\t" << gene << ":" << transcript.size() << "\n";
    }
    auto transcript_to_exon = gene_anno.transcript_to_exon;
    std::cout << "transcript_to_exon is currently " << transcript_to_exon.size() << " long\n";

    // get the junctions using transcript_to_exon
    std::unordered_map<std::string, Junctions>
    transcript_to_junctions;
    for (const auto & [transcript, exon] : transcript_to_exon) {
        transcript_to_junctions[transcript] = blocks_to_junctions(exon);
    }
    std::cout << "transcript_to_junctions is currently " << transcript_to_junctions.size() << " long\n";
    for (const auto & [transcript, junctions] : transcript_to_junctions) {
        std::cout << "\t" << transcript << ": {'right': [";
        for (const auto & junc : junctions.right) {
            std::cout << junc << ", ";
        }
        std::cout << "], 'junctions': [";
        for (const auto & junc : junctions.junctions) {
            std::cout << junc << ", ";
        }
        std::cout << "], 'left': [";
        for (const auto & junc : junctions.left) {
            std::cout << junc << ", ";
        }
        std::cout << "]}\n";
    }
    
    // remove transcripts that are too similar
    remove_similar_tr(&gene_to_transcript, &transcript_to_exon, 10);
    
    std::cout << "finished remove_rimilar_tr\n";
    auto
    gene_dict = get_gene_flat(&gene_to_transcript, &transcript_to_exon);
    std::cout << "gene_dict is currently " << gene_dict.size() << " long\n";
    
    std::cout << "finished get_gene_flat\n";
    auto
    chr_to_blocks = get_gene_blocks(&gene_dict, &chr_to_gene, &gene_to_transcript);
    std::cout << "chr_to_blocks is currently " << chr_to_blocks.size() << " long\n";
    for (const auto & [chr, block] : chr_to_blocks) {
        std::cout << "\t" << chr << ":" << block.size() << "\n";
    }

    std::cout << "finished get_gene_blocks\n";

    Config config;
    config.from_R(config_list);
    
    std::cout << "extracted config\n";
    // next we find isoforms
    std::cout << "#### Finding Isoforms\n";
    log_params(&chr_to_blocks, &gene_dict, &transcript_to_junctions, &transcript_dict);
    group_bam2isoform(
        genome_bam,
        isoform_gff3,
        tss_test_stat,
        &chr_to_blocks,
        &gene_dict,
        &transcript_to_junctions,
        &transcript_dict,
        genomefa,
        config.isoform_parameters,
        ""
    );
    std::cout << "group_bam2isoform finished\n";
    
    // get fasta
    std::cout << "doing a system command\n";
    std::string
    command = "cat " + isoform_gff3 + " > file2.gff3";
    system(command.c_str());

    std::cout << "started parsing isoform_gff3\n";
    GFFData isoform_gff = parse_gff_or_gtf(isoform_gff3);
    log_gff_data(isoform_gff, "final_gff_log_cpp.txt");

    auto chr_to_gene_iso        = isoform_gff.chr_to_gene;
    auto transcript_dict_iso    = isoform_gff.transcript_dict;
    auto gene_to_transcript_iso = isoform_gff.gene_to_transcript;
    auto transcript_to_exon_iso = isoform_gff.transcript_to_exon;

    std::cout << "assigned everything from isoform_gff\n";

    ReferenceDict
    ref_dict = {
        chr_to_gene,
        transcript_dict,
        gene_to_transcript,
        transcript_to_exon
    };

    std::cout << "wrapped them in a ReferenceDict\n";

    if (!config.realign_parameters.use_annotation) {
        ref_dict = {};
    }

    std::cout << "up to get_transcript_seq\n";
    get_transcript_seq(
        genomefa,
        transcript_fa,
        &chr_to_gene_iso,
        &transcript_dict_iso,
        &gene_to_transcript_iso,
        &transcript_to_exon_iso,
        &ref_dict
    );
    std::cout << "finished get_transcript_seq\n";

    IsoformObjects isoform_objects = {transcript_dict, transcript_dict_iso};
    
    std::ofstream
    out_test ("isoform_objects.txt");

    out_test << "transcript_dict:\n";
    for (const auto & [key, val] : isoform_objects.transcript_dict) {
        out_test << "\t" << key << ":" << val.chr << "," << val.start << "," << val.end << "\n";
    }

    out_test << "transcript_dict_iso:\n";
    for (const auto & [key, val] : isoform_objects.transcript_dict_iso) {
        out_test << "\t" << key << ":" << val.chr << "," << val.start << "," << val.end << "\n";
    }

    Rcpp::List isoform_objects_list = isoform_objects_to_R(&isoform_objects);
    return isoform_objects_to_R(&isoform_objects);
}