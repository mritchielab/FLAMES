#include "find_isoform.h"

#include "../utility/misc.h"
#include "../classes/GeneAnnoParser/GeneAnnoParser.h"
#include "../file-handling/gff3_to_fa.h"

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
    transcript_dict,

    std::string
    filename
)
{
    std::ofstream
    outfile (filename);

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
        outfile << "\t" << tr << ":{'right': " << junctions.right << ", 'junctions': (";
        if (junctions.junctions.size() > 0) {
            for (auto it = junctions.junctions.begin(); it != junctions.junctions.end() - 1; ++it) {
                outfile << *it << ", ";
            }
            outfile << junctions.junctions.back();
        }
        outfile << "), 'left': " << junctions.left << "}\n";
    }

    outfile << "\ntranscript_dict:\n";
    for (const auto & [tr, pos] : *transcript_dict) {
        outfile << "\t" << tr << ":\n\t\t(" << pos.start << "," << pos.end << ")\n";
    }
}

// [[Rcpp::export]]
Rcpp::List
find_isoform
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
    GFFData gene_anno = parseGeneAnno(gff3);
    // gene_anno.log("outfile.txt");

    // gene_anno.log("initial_gff_log_cpp_new.txt");

    auto chr_to_gene        = gene_anno.chr_to_gene;
    auto transcript_dict    = gene_anno.transcript_dict;
    auto gene_to_transcript = gene_anno.gene_to_transcript;
    auto transcript_to_exon = gene_anno.transcript_to_exon;

    // get the junctions using transcript_to_exon
    std::unordered_map<std::string, Junctions> transcript_to_junctions;
    for (const auto & [transcript, exon] : transcript_to_exon) {
        transcript_to_junctions[transcript] = blocks_to_junctions(exon); // junctions.h
    }

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
    
    std::unordered_map<std::string, std::vector<StartEndPair>>
    gene_dict = get_gene_flat(&gene_to_transcript, &transcript_to_exon); // junctions.h
    
    std::unordered_map<std::string, std::vector<GeneBlocks>>
    chr_to_blocks = get_gene_blocks(&gene_dict, &chr_to_gene, &gene_to_transcript); // junctions.h
    
    Config config;
    config.from_R(config_list);
    
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

    GFFData isoform_gff = parseGeneAnno(isoform_gff3);
    // isoform_gff.log("final_gff_log_cpp.txt");

    auto chr_to_gene_iso        = isoform_gff.chr_to_gene;
    auto transcript_dict_iso    = isoform_gff.transcript_dict;
    auto gene_to_transcript_iso = isoform_gff.gene_to_transcript;
    auto transcript_to_exon_iso = isoform_gff.transcript_to_exon;

    ReferenceDict
    ref_dict = {
        chr_to_gene,
        transcript_dict,
        gene_to_transcript,
        transcript_to_exon
    };

    if (!config.realign_parameters.use_annotation) {
        ref_dict = {};
    }

    get_transcript_seq(
        genomefa,
        transcript_fa,
        &chr_to_gene_iso,
        &transcript_dict_iso,
        &gene_to_transcript_iso,
        &transcript_to_exon_iso,
        &ref_dict
    );

    IsoformObjects isoform_objects = {transcript_dict, transcript_dict_iso};
    
    // std::ofstream
    // out_test ("isoform_objects.txt");

    // out_test << "transcript_dict:\n";
    // for (const auto & [key, val] : isoform_objects.transcript_dict) {
    //     out_test << "\t" << key << ":" << val.chr << "," << val.start << "," << val.end << "\n";
    // }

    // out_test << "transcript_dict_iso:\n";
    
    Rcpp::List isoform_objects_list = isoform_objects_to_R(&isoform_objects);
    return isoform_objects_to_R(&isoform_objects);
}