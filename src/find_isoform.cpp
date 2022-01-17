#include <string>
#include <unordered_map>
#include <vector>
#include <Rcpp.h>

#include "config.h"
#include "misc.h"
#include "parse_gene_anno_native.h"
#include "junctions.h"
#include "gff3_to_fa.hpp"
#include "group_bam2isoform.h"
#include "ReferenceDict.hpp"

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

// [[Rcpp::export]]
Rcpp::List
find_isoform_cpp
(
    std::string gff3, 
    std::string genome_bam, 
    std::string isoform_gff3, 
    std::string tss_tes_stat, 
    std::string genomefa, 
    std::string transcript_fa, 
    int         downsample_ratio, 
    Rcpp::List  config_list, 
    std::string raw_splice_isoform
)
{
    Rcpp::Rcout << "#### Reading Gene Annotations\n";

    // first, extract from the gff file
    GFFData gene_anno = parse_gff_or_gtf(gff3); // parse_gene_anno_native.h

    std::unordered_map<std::string, std::vector<std::string>> chr_to_gene        = gene_anno.chr_to_gene;
    Rcpp::Rcout << "chr_to_gene is size " << chr_to_gene.size() << "\n";
    for (const auto & [chr, gene] : chr_to_gene) {
        Rcpp::Rcout << "\t" << chr << ":" << gene.size() << "\n";
    }
    Rcpp::Rcout << "chr_to_gene is currently " << chr_to_gene.size() << " long\n";

    std::unordered_map<std::string, Pos> transcript_dict    = gene_anno.transcript_dict;
    Rcpp::Rcout << "transcript_dict is currently " << transcript_dict.size() << " long\n";

    std::unordered_map<std::string, std::vector<std::string>> gene_to_transcript = gene_anno.gene_to_transcript;
    Rcpp::Rcout << "gene_to_transcript is currently " << gene_to_transcript.size() << " long\n";
    for (const auto & [gene, transcript] : gene_to_transcript) {
        Rcpp::Rcout << "\t" << gene << ":" << transcript.size() << "\n";
    }
    std::unordered_map<std::string, std::vector<StartEndPair>> transcript_to_exon = gene_anno.transcript_to_exon;
    Rcpp::Rcout << "transcript_to_exon is currently " << transcript_to_exon.size() << " long\n";

    // get the junctions using transcript_to_exon
    std::unordered_map<std::string, Junctions> transcript_to_junctions;
    for (const auto & [transcript, exon] : transcript_to_exon) {
        transcript_to_junctions[transcript] = blocks_to_junctions(exon); // junctions.h
    }
    Rcpp::Rcout << "transcript_to_junctions is currently " << transcript_to_junctions.size() << " long\n";


    // Rcpp::Rcout << "gene_to_transcript:\n";
    // for (const auto & [gene, transcript] : gene_to_transcript) {
    //     Rcpp::Rcout << "gene:" << gene << "\ntranscript:";
    //     for (const auto & tr : transcript) {
    //         Rcpp::Rcout << tr << ",";
    //     }
    // }

    // Rcpp::Rcout << "transcript_to_junctions:\n";
    // for (const auto & [transcript, junctions] : transcript_to_junctions) {
    //     Rcpp::Rcout << "transcript:" << transcript << "\n";
    //     Rcpp::Rcout << "junctions:" << junctions.left[0] << "," << junctions.right[0] << "\n";
    // }
    
    // remove transcripts that are too similar
    remove_similar_tr(gene_to_transcript, transcript_to_exon, 10); // junctions.h
    
    Rcpp::Rcout << "finished remove_rimilar_tr\n";
    std::unordered_map<std::string, std::vector<StartEndPair>>
    gene_dict = get_gene_flat(&gene_to_transcript, &transcript_to_exon); // junctions.h
    Rcpp::Rcout << "gene_dict is currently " << gene_dict.size() << " long\n";
    
    Rcpp::Rcout << "finished get_gene_flat\n";
    std::unordered_map<std::string, std::vector<GeneBlocks>>
    chr_to_blocks = get_gene_blocks(&gene_dict, &chr_to_gene, &gene_to_transcript); // junctions.h
    Rcpp::Rcout << "chr_to_blocks is currently " << chr_to_blocks.size() << " long\n";
    for (const auto & [chr, block] : chr_to_blocks) {
        Rcpp::Rcout << "\t" << chr << ":" << block.size() << "\n";
    }

    Rcpp::Rcout << "finished get_gene_blocks\n";

    Config config;
    config.from_R(config_list);
    
    Rcpp::Rcout << "extracted config\n";
    // next we find isoforms
    Rcpp::Rcout << "#### Finding Isoforms\n";
    group_bam2isoform(
    // minimal_group_bam2isoform(
        genome_bam,
        isoform_gff3,
        tss_tes_stat,
        &chr_to_blocks,
        &gene_dict,
        &transcript_to_junctions,
        &transcript_dict,
        genomefa,
        config.isoform_parameters,
        ""
    );

    Rcpp::Rcout << "group_bam2isoform finished\n";

    // get fasta
    GFFData isoform_gff = parse_gff_or_gtf(isoform_gff3);

    auto chr_to_gene_iso        = isoform_gff.chr_to_gene;
    auto transcript_dict_iso    = isoform_gff.transcript_dict;
    auto gene_to_transcript_iso = isoform_gff.gene_to_transcript;
    auto transcript_to_exon_iso = isoform_gff.transcript_to_exon;

    Rcpp::Rcout << "assigned everything from isoform_gff\n";

    ReferenceDict
    ref_dict = {
        chr_to_gene,
        transcript_dict,
        gene_to_transcript,
        transcript_to_exon
    };

    Rcpp::Rcout << "wrapped them in a ReferenceDict\n";

    if (!config.realign_parameters.use_annotation) {
        ref_dict = {};
    }

    Rcpp::Rcout << "up to get_transcript_seq\n";
    get_transcript_seq(
        genomefa,
        transcript_fa,
        &chr_to_gene_iso,
        &transcript_dict_iso,
        &gene_to_transcript_iso,
        &transcript_to_exon_iso,
        &ref_dict
    );
    Rcpp::Rcout << "finished get_transcript_seq\n";

    IsoformObjects isoform_objects = {transcript_dict, transcript_dict_iso};
    
    std::ofstream
    out_test ("isoform_objects.txt");

    out_test << "transcript_dict:\n";
    for (const auto & [key, val] : isoform_objects.transcript_dict) {
        out_test << "\t" << key << ":" << val.chr << "," << val.start << "," << val.end << "\n";
    }

    out_test << "transcript_dict_iso:\n";
    

    Rcpp::List isoform_objects_list = isoform_objects_to_R(&isoform_objects);
    // std::ofstream
    // out_test2 ("isoform_objects_r.txt");
    // for (const auto & [key, val] : isoform_objects_list["transcript_dict"]) {
    //     out_test << "\t" << key << ":" << val.chr << "," << val.start << "," << val.end << "\n";
    // }



    return isoform_objects_to_R(&isoform_objects);
}