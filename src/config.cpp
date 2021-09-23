#include "config.h"
#include <Rcpp.h>

namespace Rcpp {
    template <>
    SEXP wrap(const PipelineParameters& pipeline_parameters)
    {
        return Rcpp::wrap(Rcpp::List::create(
            Rcpp::Named("do_genome_alignment") = Rcpp::wrap(pipeline_parameters.do_genome_alignment),
            Rcpp::Named("do_isoform_identification") = Rcpp::wrap(pipeline_parameters.do_isoform_identification),
            Rcpp::Named("do_read_realignment") = Rcpp::wrap(pipeline_parameters.do_read_realignment),
            Rcpp::Named("do_transcript_quantification") = Rcpp::wrap(pipeline_parameters.do_transcript_quantification)
        ));
    };

    template <>
    SEXP wrap(const GlobalParameters& global_parameters)
    {
        return Rcpp::wrap(Rcpp::List::create(
            Rcpp::Named("generate_raw_isoform") = Rcpp::wrap(global_parameters.generate_raw_isoform),
            Rcpp::Named("has_UMI") = Rcpp::wrap(global_parameters.has_UMI)
        ));
    };

    template <>
    SEXP wrap(const IsoformParameters& isoform_parameters)
    {
        return Rcpp::wrap(Rcpp::List::create(
            Rcpp::Named("MAX_DIST") = Rcpp::wrap(isoform_parameters.MAX_DIST),
            Rcpp::Named("MAX_TS_DIST") = Rcpp::wrap(isoform_parameters.MAX_TS_DIST),
            Rcpp::Named("MAX_SPLICE_MATCH_DIST") = Rcpp::wrap(isoform_parameters.MAX_SPLICE_MATCH_DIST),
            Rcpp::Named("MIN_FL_EXON_LEN") = Rcpp::wrap(isoform_parameters.MIN_FL_EXON_LEN),
            Rcpp::Named("MAX_SITE_PER_SPLICE") = Rcpp::wrap(isoform_parameters.MAX_SITE_PER_SPLICE),
            Rcpp::Named("MIN_SUP_CNT") = Rcpp::wrap(isoform_parameters.MIN_SUP_CNT),
            
            Rcpp::Named("MIN_CNT_PCT") = Rcpp::wrap(isoform_parameters.MIN_CNT_PCT),
            Rcpp::Named("MIN_SUP_PCT") = Rcpp::wrap(isoform_parameters.MIN_SUP_PCT),

            Rcpp::Named("STRAND_SPECIFIC") = Rcpp::wrap(isoform_parameters.STRAND_SPECIFIC),
            Rcpp::Named("REMOVE_INCOMP_READS") = Rcpp::wrap(isoform_parameters.REMOVE_INCOMP_READS)
        ));
    };

    template <>
    SEXP wrap(const AlignmentParameters& alignment_parameters) {
        return Rcpp::wrap(Rcpp::List::create(
            Rcpp::Named("use_junctions") = Rcpp::wrap(alignment_parameters.use_junctions),
            Rcpp::Named("no_flank") = Rcpp::wrap(alignment_parameters.no_flank)
        ));
    };

    template <>
    SEXP wrap(const RealignParameters& realign_parameters) {
        return Rcpp::wrap(Rcpp::List::create(
            Rcpp::Named("use_annotation") = Rcpp::wrap(realign_parameters.use_annotation)
        ));
    };
    
    template <>
    SEXP wrap(const TranscriptCounting& transcript_counting) {
        return Rcpp::wrap(Rcpp::List::create(
            Rcpp::Named("min_tr_coverage") = Rcpp::wrap(transcript_counting.min_tr_coverage),
            Rcpp::Named("min_read_coverage") = Rcpp::wrap(transcript_counting.min_read_coverage)
        ));
    };
    
    template <>
    SEXP wrap(const Config& config) {
        return Rcpp::wrap(Rcpp::List::create(
            Rcpp::Named("pipeline_parameters") = Rcpp::wrap(config.pipeline_parameters),
            Rcpp::Named("global_parameters") = Rcpp::wrap(config.global_parameters),
            Rcpp::Named("isoform_parameters") = Rcpp::wrap(config.isoform_parameters),
            Rcpp::Named("alignment_parameters") = Rcpp::wrap(config.alignment_parameters),
            Rcpp::Named("realign_parameters") = Rcpp::wrap(config.realign_parameters),
            Rcpp::Named("transcript_counting") = Rcpp::wrap(config.transcript_counting)
        ));
    };
}