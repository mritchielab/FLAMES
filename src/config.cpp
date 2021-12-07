#include <Rcpp.h>

#include "config.h"

Rcpp::List
PipelineParameters::to_R()
{
    return Rcpp::wrap(Rcpp::List::create(
        Rcpp::Named("do_genome_alignment") = Rcpp::wrap(this->do_genome_alignment),
        Rcpp::Named("do_isoform_identification") = Rcpp::wrap(this->do_isoform_identification),
        Rcpp::Named("do_read_realignment") = Rcpp::wrap(this->do_read_realignment),
        Rcpp::Named("do_transcript_quantification") = Rcpp::wrap(this->do_transcript_quantification)
    ));
}

void
PipelineParameters::from_R(Rcpp::List list)
{
    this->do_genome_alignment = list["do_genome_alignment"];
    this->do_isoform_identification = list["do_isoform_identification"];
    this->do_read_realignment = list["do_read_realignment"];
    this->do_transcript_quantification = list["do_transcript_quantification"];
}

void
PipelineParameters::print()
{
  Rcpp::Rcout << "\tpipeline_parameters\n"
            << "\t\tdo_genome_alignment : " << this->do_genome_alignment << "\n"
            << "\t\tdo_isoform_identification : " << this->do_isoform_identification << "\n"
            << "\t\tdo_read_realignment : " << this->do_read_realignment << "\n"
            << "\t\tdo_transcript_quantification : " << this->do_transcript_quantification << "\n";
}

Rcpp::List
GlobalParameters::to_R()
{
    return Rcpp::wrap(Rcpp::List::create(
        Rcpp::Named("generate_raw_isoform") = Rcpp::wrap(this->generate_raw_isoform),
        Rcpp::Named("has_UMI") = Rcpp::wrap(this->has_UMI)
    ));
}
 
void
GlobalParameters::from_R(Rcpp::List list)
{
    this->generate_raw_isoform = list["generate_raw_isoform"];
    this->has_UMI = list["has_UMI"];
}

void
GlobalParameters::print()
{
  Rcpp::Rcout << "\tglobal parameters\n"
            << "\t\tgenerate_raw_isoform : " << this->generate_raw_isoform << "\n"
            << "\t\thas_UMI : " << this->has_UMI << "\n";
}

Rcpp::List
IsoformParameters::to_R()
{
    return Rcpp::wrap(Rcpp::List::create(
        Rcpp::Named("MAX_DIST") = Rcpp::wrap(this->MAX_DIST),
        Rcpp::Named("MAX_TS_DIST") = Rcpp::wrap(this->MAX_TS_DIST),
        Rcpp::Named("MAX_SPLICE_MATCH_DIST") = Rcpp::wrap(this->MAX_SPLICE_MATCH_DIST),
        Rcpp::Named("MIN_FL_EXON_LEN") = Rcpp::wrap(this->MIN_FL_EXON_LEN),
        Rcpp::Named("MAX_SITE_PER_SPLICE") = Rcpp::wrap(this->MAX_SITE_PER_SPLICE),
        Rcpp::Named("MIN_SUP_CNT") = Rcpp::wrap(this->MIN_SUP_CNT),
        
        Rcpp::Named("MIN_CNT_PCT") = Rcpp::wrap(this->MIN_CNT_PCT),
        Rcpp::Named("MIN_SUP_PCT") = Rcpp::wrap(this->MIN_SUP_PCT),

        Rcpp::Named("STRAND_SPECIFIC") = Rcpp::wrap(this->STRAND_SPECIFIC),
        Rcpp::Named("REMOVE_INCOMP_READS") = Rcpp::wrap(this->REMOVE_INCOMP_READS)
    ));
}

void
IsoformParameters::from_R(Rcpp::List list)
{
    this->MAX_DIST                = list["MAX_DIST"];
    this->MAX_TS_DIST             = list["MAX_TS_DIST"];
    this->MAX_SPLICE_MATCH_DIST   = list["MAX_SPLICE_MATCH_DIST"];
    this->MIN_FL_EXON_LEN         = list["MIN_FL_EXON_LEN"];
    this->MAX_SITE_PER_SPLICE     = list["MAX_SITE_PER_SPLICE"];
    this->MIN_SUP_CNT             = list["MIN_SUP_CNT"];
    
    this->MIN_CNT_PCT             = list["MIN_CNT_PCT"];
    this->MIN_SUP_PCT             = list["MIN_SUP_PCT"];

    this->STRAND_SPECIFIC         = list["STRAND_SPECIFIC"];
    this->REMOVE_INCOMP_READS     = list["REMOVE_INCOMP_READS"];
}

void
IsoformParameters::print()
{
  Rcpp::Rcout << "\tisoform_parameters\n"
            << "\t\tMAX_DIST : " << this->MAX_DIST << "\n"
            << "\t\tMAX_TS_DIST : " << this->MAX_TS_DIST << "\n"
            << "\t\tMAX_SPLICE_MATCH_DIST : " << this->MAX_SPLICE_MATCH_DIST << "\n"
            << "\t\tMIN_FL_EXON_LEN : " << this->MIN_FL_EXON_LEN << "\n"
            << "\t\tMAX_SITE_PER_SPLICE : " << this->MAX_SITE_PER_SPLICE << "\n"
            << "\t\tMIN_SUP_CNT : " << this->MIN_SUP_CNT << "\n"

            << "\t\tMIN_CNT_PCT : " << this->MIN_CNT_PCT << "\n"
            << "\t\tMIN_SUP_PCT : " << this->MIN_SUP_PCT << "\n"

            << "\t\tSTRAND_SPECIFIC : " << this->STRAND_SPECIFIC << "\n"
            << "\t\tREMOVE_INCOMP_READS : " << this->REMOVE_INCOMP_READS << "\n";
}

Rcpp::List
AlignmentParameters::to_R()
{
    return Rcpp::wrap(Rcpp::List::create(
         Rcpp::Named("use_junctions") = Rcpp::wrap(this->use_junctions),
         Rcpp::Named("no_flank") = Rcpp::wrap(this->no_flank)
    ));
}

void
AlignmentParameters::from_R(Rcpp::List list)
{
    this->use_junctions = list["use_junctions"];
    this->no_flank = list["no_flank"];
}

void
AlignmentParameters::print()
{
    Rcpp::Rcout << "\talignment parameters\n"
              << "\t\tuse_junctions : " << this->use_junctions << "\n"
              << "\t\tno_flank : " << this->no_flank << "\n";
}

Rcpp::List
RealignParameters::to_R()
{
    return Rcpp::wrap(Rcpp::List::create(
         Rcpp::Named("use_annotation") = Rcpp::wrap(this->use_annotation)
    ));
}

void
RealignParameters::from_R(Rcpp::List list)
{
    this->use_annotation = list["use_annotation"];
}

void
RealignParameters::print()
{
    Rcpp::Rcout << "\trealign parameters\n"
              << "\t\tuse_annotation : " << this->use_annotation << "\n";
}

Rcpp::List
TranscriptCounting::to_R()
{
    return Rcpp::wrap(Rcpp::List::create(
        Rcpp::Named("min_tr_coverage") = Rcpp::wrap(this->min_tr_coverage),
        Rcpp::Named("min_read_coverage") = Rcpp::wrap(this->min_read_coverage)
    ));
}

void
TranscriptCounting::from_R(Rcpp::List list)
{
    this->min_tr_coverage = list["min_tr_coverage"];
    this->min_read_coverage = list["min_read_coverage"];
}

void
TranscriptCounting::print()
{
    Rcpp::Rcout << "\ttranscript counting\n"
              << "\t\tmin_tr_coverage : " << this->min_tr_coverage << "\n"
              << "\t\tmin_read_coverage : " << this->min_read_coverage << "\n";
}

Rcpp::List
Config::to_R()
{
    return Rcpp::wrap(Rcpp::List::create(
        Rcpp::Named("pipeline_parameters") = Rcpp::wrap(this->pipeline_parameters.to_R()),
        Rcpp::Named("global_parameters") = Rcpp::wrap(this->global_parameters.to_R()),
        Rcpp::Named("isoform_parameters") = Rcpp::wrap(this->isoform_parameters.to_R()),
        Rcpp::Named("alignment_parameters") = Rcpp::wrap(this->alignment_parameters.to_R()),
        Rcpp::Named("realign_parameters") = Rcpp::wrap(this->realign_parameters.to_R()),
        Rcpp::Named("transcript_counting") = Rcpp::wrap(this->transcript_counting.to_R())
    ));
}

void
Config::from_R(Rcpp::List list)
{
    this->pipeline_parameters.from_R(list["pipeline_parameters"]);
    this->global_parameters.from_R(list["global_parameters"]);
    this->alignment_parameters.from_R(list["alignment_parameters"]);
    this->realign_parameters.from_R(list["realign_parameters"]);
    this->transcript_counting.from_R(list["transcript_counting"]);
}

Config::Config(Rcpp::List list) {
	from_R(list);
}
Config::Config() {
}

void
Config::print()
{
    this->pipeline_parameters.print();
    this->global_parameters.print();
    this->isoform_parameters.print();
    this->alignment_parameters.print();
    this->realign_parameters.print();
    this->transcript_counting.print();
}
