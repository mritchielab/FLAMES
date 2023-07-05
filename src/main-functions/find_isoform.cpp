#include "find_isoform.h"

#include <string>
#include <unordered_map>
#include <vector>
#include <utility>

#include "../classes/GFFData.h"
#include "../classes/GeneAnnotationParser.h"
#include "../classes/junctions.h"
#include "../classes/Pos.h"
#include "get_transcript_seq.h"
#include "group_bam2isoform.h"

void
find_isoform
(
    const std::string &gff3,
    const std::string &genome_bam,
    const std::string &isoform_gff3,
    const std::string &tss_tes_stat,
    const std::string &genomefa,
    const std::string &transcript_fa,
    const Rcpp::List  &isoform_parameters,
    const std::string &raw_splice_isoform)
{
    Rcpp::Rcout << "#### Reading Gene Annotations\n";

    GFFData gene_annotation = parse_gff_file(gff3);

    std::unordered_map<std::string, Junctions>
    transcript_to_junctions = map_transcripts_to_junctions(
        gene_annotation.transcript_to_exon
    );

    gene_annotation.gene_to_transcript = remove_similar_tr(
        gene_annotation.gene_to_transcript,
        gene_annotation.transcript_to_exon,
        10
    );

    std::unordered_map<std::string, std::vector<exon>>
    gene_dict = get_gene_flat(
        gene_annotation.gene_to_transcript,
        gene_annotation.transcript_to_exon
    );

    std::unordered_map<std::string, std::vector<GeneBlocks>>
    chr_to_blocks = get_gene_blocks(
        gene_dict, 
        gene_annotation.chr_to_gene,
        gene_annotation.gene_to_transcript
    );

    // GROUP_BAM2ISOFORM
    group_bam2isoform(
        genome_bam,
        isoform_gff3,
        tss_tes_stat,
        gene_dict,
        transcript_to_junctions,
        gene_annotation.transcript_dict,
        chr_to_blocks,
        genomefa,
        isoform_parameters,
        raw_splice_isoform
    );

    GFFData isoform_annotation = parse_gff_file(isoform_gff3);
    
    // get_transcript_seq
    get_transcript_seq(
        genomefa,
        transcript_fa,
        isoform_annotation,
        gene_annotation
    ); // - This does not modify values that are used later (it modifies chr_to_blocks, but only transcript_dict_i and transcript_dict are used later.)
    return;
}