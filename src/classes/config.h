#ifndef CONFIG_H
#define CONFIG_H

#include <Rcpp.h>

class PipelineParameters
{
    public:
        bool do_genome_alignment = true;
        bool do_isoform_identification = true;
        bool do_read_realignment = true;
        bool do_transcript_quantification = true;

        Rcpp::List 
        to_R();
        void 
        from_R(Rcpp::List list); 
        void 
        print();
};

class GlobalParameters
{
    public:
        bool generate_raw_isoform = true;
        bool has_UMI = false;

        Rcpp::List to_R();
        void from_R(Rcpp::List list); 
        void print();
};

class IsoformParameters 
{
    public:
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

        Rcpp::List to_R();
        void from_R(Rcpp::List list); 
        void print();
};

class AlignmentParameters 
{
    public:
        bool use_junctions = true;
        bool no_flank = true;

        Rcpp::List to_R();
        void from_R(Rcpp::List list); 
        void print();
};

class RealignParameters
{
    public:
        bool use_annotation = true;

        Rcpp::List to_R();
        void from_R(Rcpp::List list); 
        void print();
};

class TranscriptCounting
{
    public:
        float min_tr_coverage = 0.75;
        float min_read_coverage = 0.75;

        Rcpp::List to_R();
        void from_R(Rcpp::List list); 
        void print();
};

class Config
{
    public:
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

        Rcpp::List 
        to_R();

        void
        from_R(Rcpp::List list);

		Config();
		
        Config(Rcpp::List list);
        
        void
        print();
};

#endif // CONFIG_H