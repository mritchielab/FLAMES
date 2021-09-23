#include <string>
#include <vector>
#include <map>
#include <RcppCommon.h>

#ifndef CONFIG
#define CONFIG

struct PipelineParameters
{
    bool do_genome_alignment = true;
    bool do_isoform_identification = true;
    bool do_read_realignment = true;
    bool do_transcript_quantification = true;
};

struct GlobalParameters
{
    bool generate_raw_isoform = true;
    bool has_UMI = false;
};

struct IsoformParameters 
{
    int MAX_DIST = 10;
    int MAX_TS_DIST = 100;
    int MAX_SPLICE_MATCH_DIST = 10;
    int MIN_FL_EXON_LEN = 40;
    int MAX_SITE_PER_SPLICE = 3;
    int MIN_SUP_CNT = 10;

    float MIN_CNT_PCT = 0.01;
    float MIN_SUP_PCT = 0.2;

    int STRAND_SPECIFIC = 1;
    int REMOVE_INCOMP_READS = 5;
};

struct AlignmentParameters 
{
    bool use_junctions = true;
    bool no_flank = true;
};

struct RealignParameters
{
    bool use_annotation = true;
};

struct TranscriptCounting
{
    float min_tr_coverage = 0.75;
    float min_read_coverage = 0.75;
};

struct Config
{
    PipelineParameters
    pipeline_parameters;

    GlobalParameters
    global_parameters;

    IsoformParameters
    isoform_parameters;

    AlignmentParameters
    alignment_parameters;

    RealignParameters
    realign_parameters;

    TranscriptCounting
    transcript_counting;
};

#endif

namespace Rcpp {
    template <>
    SEXP wrap(const PipelineParameters& pipeline_parameters);

    template <>
    SEXP wrap(const GlobalParameters& global_parameters);

    template <>
    SEXP wrap(const IsoformParameters& isoform_parameters);

    template <>
    SEXP wrap(const AlignmentParameters& alignment_parameters);

    template <>
    SEXP wrap(const RealignParameters& realign_parameters);
    
    template <>
    SEXP wrap(const TranscriptCounting& transcript_counting);
    
    template <>
    SEXP wrap(const Config& config);
}