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

    auto chr_to_gene        = gene_anno.chr_to_gene;
    auto transcript_dict    = gene_anno.transcript_dict;
    auto gene_to_transcript = gene_anno.gene_to_transcript;
    auto transcript_to_exon = gene_anno.transcript_to_exon;

    // get the junctions using transcript_to_exon
    std::unordered_map<std::string, Junctions>
    transcript_to_junctions;
    for (const auto & [transcript, exon] : transcript_to_exon) {
        transcript_to_junctions[transcript] = blocks_to_junctions(exon);
    }

    // std::cout << "gene_to_transcript:\n";
    // for (const auto & [gene, transcript] : gene_to_transcript) {
    //     std::cout << "gene:" << gene << "\ntranscript:";
    //     for (const auto & tr : transcript) {
    //         std::cout << tr << ",";
    //     }
    // }

    // std::cout << "transcript_to_junctions:\n";
    // for (const auto & [transcript, junctions] : transcript_to_junctions) {
    //     std::cout << "transcript:" << transcript << "\n";
    //     std::cout << "junctions:" << junctions.left[0] << "," << junctions.right[0] << "\n";
    // }
    
    // remove transcripts that are too similar
    remove_similar_tr(&gene_to_transcript, &transcript_to_exon, 10);
    
    std::cout << "finished remove_rimilar_tr\n";
    auto
    gene_dict = get_gene_flat(&gene_to_transcript, &transcript_to_exon);
    
    std::cout << "finished get_gene_flat\n";
    auto
    chr_to_blocks = get_gene_blocks(&gene_dict, &chr_to_gene, &gene_to_transcript);

    std::cout << "finished get_gene_blocks\n";

    Config config;
    config.from_R(config_list);
    
    std::cout << "extracted config\n";
    // next we find isoforms
    std::cout << "#### Finding Isoforms\n";
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
    GFFData isoform_gff = parse_gff_or_gtf(isoform_gff3);

    auto chr_to_gene_iso        = isoform_gff.chr_to_gene;
    auto transcript_dict_iso    = isoform_gff.transcript_dict;
    auto gene_to_transcript_iso = isoform_gff.gene_to_transcript;
    auto transcript_to_exon_iso = isoform_gff.transcript_to_exon;

    std::cout << "assigned everything from isoform_gff\n";
    
    std::cout << "chr_to_gene:\n";
    for (const auto & [chr, gene] : chr_to_gene) {
        std::cout << "\tchr:" << chr << "\n";
        std::cout << "\tgene:" << gene.front() << " ... " << gene.back() << "\n";
    }

    std::cout << "transcript_dict:\n";
    for (const auto & [tr, pos] : transcript_dict) {
        std::cout << "\ttr:" << tr << "\n";
        std::cout << "\tpos:" << pos.start << ", " << pos.end << "\n";
    }

    std::cout << "gene_to_transcript:\n";
    for (const auto & [gene, transcript] : gene_to_transcript) {
        std::cout << "\tgene:" << gene << "\n";
        std::cout << "\ttranscript:" << transcript.front() << "\n";
    }

    std::cout << "transcript_to_exon:\n";
    for (const auto & [tr, exon] : transcript_to_exon) {
        std::cout << "\ttr:" << tr << "\n";
        std::cout << "\tex:" << exon.front().start << "\n";
    }

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

    return isoform_objects_to_R(&isoform_objects);
}