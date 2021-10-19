#include "find_isoform.h"

// [[Rcpp::export]]
void
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
    std::cout << "### Reading Gene Annotations\n";

    // first, extract from the gff file
    GFFData gene_anno = parse_gff_tree(gff3);

    auto chr_to_gene        = gene_anno.chr_to_gene;
    auto transcript_dict    = gene_anno.transcript_dict;
    auto gene_to_transcript = gene_anno.gene_to_transcript;
    auto transcript_to_exon = gene_anno.transcript_to_exon;

    // convert the chr_to_gene list to a std::unordered_map
    // convert the gene_to_transcript list to a std::unordered_map

    // get the junctions using transcript_to_exon

    std::unordered_map<std::string, Junctions>
    transcript_to_junctions;
    for (const auto & [transcript, exon] : transcript_to_exon) {
        transcript_to_junctions[transcript] = blocks_to_junctions(exon);
    }

    // remove transcripts that are too similar
    remove_similar_tr(&gene_to_transcript, &transcript_to_exon, 10);
    
    auto
    gene_dict = get_gene_flat(&gene_to_transcript, &transcript_to_exon);

    auto
    chr_to_blocks = get_gene_blocks(&gene_dict, &chr_to_gene, &gene_to_transcript);

    Config config;
    config.from_R(config_list);

    // next we find isoforms
    std::cout << "### Finding Isoforms\n";
    // group_bam2isoform(
    //     genome_bam,
    //     isoform_gff3,
    //     tss_test_stat,
    //     "",
    //     &chr_to_blocks,
    //     &gene_dict,
    //     &transcript_to_junctions,
    //     &transcript_dict,
    //     &genomefa,
    //     &(config.isoform_parameters)
    // );

    // get fasta
    GFFData isoform_gff = parse_gff_tree(isoform_gff3);

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

    get_transcript_seq(
        genomefa,
        transcript_fa,
        &chr_to_gene_iso,
        &transcript_dict_iso,
        &gene_to_transcript_iso,
        &transcript_to_exon_iso,
        &ref_dict
    );
}