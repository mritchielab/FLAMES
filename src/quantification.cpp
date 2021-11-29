#include "quantification.hpp"

// [[Rcpp::export]]
void
quantification_cpp
(
    Rcpp::List  config_list,
    std::string realign_bam,
    std::string transcript_fa_idx,
    Rcpp::List  isoform_objects_list,
    std::string tr_cnt_csv,
    std::string tr_badcov_cnt_csv,
    std::string isoform_gff3,
    std::string annot,
    std::string isoform_gff3_f,
    std::string FSM_anno_out
)
{
    std::cout << "#### Generating transcript count matrix\n";
    
    GFFData gene_anno = parse_gff_or_gtf(annot);
    auto transcript_to_exon = gene_anno.transcript_to_exon;
    std::cout << "transcript_to_exon is currently " << transcript_to_exon.size() << " isoforms long\n";

    // unwrap things from R
    Config config;
    config.from_R(config_list);
    IsoformObjects isoform_objects = isoform_objects_from_R(isoform_objects_list);
    std::cout << "isoform_objects.transcript_dict is " << isoform_objects.transcript_dict.size() << " long\n";
    std::cout << "isoform_objects.transcript_dict_iso is " << isoform_objects.transcript_dict_iso.size() << " long\n";
    
    // std::cout << "isoform_objects.transcript_dict:\n";
    // for (const auto & [key, val] : isoform_objects.transcript_dict) {
    //     std::cout << "\t" << key << ":" << val.chr << "," << val.start << "," << val.end << "\n";
    // }
    // return;

    auto
    parse_realign = parse_realigned_bam(
        realign_bam,
        transcript_fa_idx,
        config.isoform_parameters.MIN_SUP_CNT,
        config.transcript_counting.min_tr_coverage,
        config.transcript_counting.min_read_coverage,
        {}
    );

    auto
    tr_count = write_tr_to_csv_cpp(
        parse_realign.bc_tr_count_dict,
        isoform_objects.transcript_dict_iso,
        tr_cnt_csv,
        isoform_objects.transcript_dict,
        config.global_parameters.has_UMI
    );

    write_tr_to_csv_cpp(
        parse_realign.bc_tr_badcov_count_dict,
        isoform_objects.transcript_dict_iso,
        tr_badcov_cnt_csv,
        isoform_objects.transcript_dict,
        config.global_parameters.has_UMI
    );

    annotate_filter_gff(
        isoform_gff3,
        annot,
        isoform_gff3_f,
        FSM_anno_out,
        tr_count,
        config.isoform_parameters.MIN_SUP_CNT
    );
}