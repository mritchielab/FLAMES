#include <string>
#include <iostream>
#include <Rcpp.h>
#include <fstream>
#include <array>
#include <stdexcept>
#include <assert.h>

#include "../utility/json/json.h"
#include "../classes/Config.h"
#include "parse_json_config.h"

// [[Rcpp::export]]
Rcpp::List
parse_json_config(std::string json_file) {
    std::ifstream file(json_file);

    Json::Value json;

    file >> json;

    // Rcpp::Rcout << json["windowed"] << "\n" << json["res"][0] << "\n";
    if (!verify_json_config(json)) {
        Rcpp::Rcout << "problem with config; expect errors!!\n";
    }

    return load_json_config(json).to_R();
}

// [[Rcpp::export]]
void print_config(Rcpp::List list) {
    Config config(list);
    Rcpp::Rcout << "\tParameters in configuration file:\n";
    config.print();
}

int verify_json_config(Json::Value json) {
    // check that the json itself is valid
    if (!(json.isObject())) {
        return 0;
    }

    // check all base parameters
    std::array < std::string, 5 > parameters = {
        "transcript_counting",
        "isoform_parameters",
        "alignment_parameters",
        "global_parameters",
        "pipeline_parameters"
    };
    for (auto parameter: parameters) {
        if (!(json.isMember(parameter))) {
            return 0;
        }
    }

    // check the contents of isoform_parameters
    std::array < std::string, 9 > isoform_parameters = {
        "MAX_DIST",
        "MAX_TS_DIST",
        "MAX_SPLICE_MATCH_DIST",
        "Max_site_per_splice",
        "Min_sup_cnt",
        "Min_cnt_pct",
        "Min_sup_pct",
        "strand_specific",
        "remove_incomp_reads"
    };
    for (auto parameter: isoform_parameters) {
        if (!(json["isoform_parameters"].isMember(parameter))) {
            return 0;
        }
    }

    // check the contents of alignment_parameters
    if (!(json["alignment_parameters"].isMember("use_junctions"))) {
        return 0;
    }

    // check the contents of alignment_parameters
    if (!(json["global_parameters"].isMember("generate_raw_isoform"))) {
        return 0;
    }

    // check the contents of pipeline_parameters
    std::array < std::string, 4 > pipeline_parameters = {
        "do_genome_alignment",
        "do_isoform_identification",
        "do_read_realignment",
        "do_transcript_quantification"
    };
    for (auto parameter: pipeline_parameters) {
        if (!(json["pipeline_parameters"].isMember(parameter))) {
            return 0;
        }
    }

    // check the contents of transcript_counting
    std::array < std::string, 2 > transcript_counting = {
        "min_tr_coverage",
        "min_read_coverage"
    };
    for (auto parameter: transcript_counting) {
        if (!(json["transcript_counting"].isMember(parameter))) {
            return 0;
        }
    }

    // check all arguments are within the range
    if (!(json["isoform_parameters"]["MAX_DIST"] > 0 &&
            json["isoform_parameters"]["MAX_TS_DIST"] > 0 &&
            json["isoform_parameters"]["MAX_SPLICE_MATCH_DIST"] > 0 &&
            json["isoform_parameters"]["Max_site_per_splice"] > 0 &&
            json["isoform_parameters"]["Min_sup_cnt"] > 0 &&
            json["isoform_parameters"]["Min_sup_pct"] > 0.0 &&
            json["isoform_parameters"]["Min_sup_pct"] < 1.0 &&
            (json["isoform_parameters"]["strand_specific"] == -1 ||
                json["isoform_parameters"]["strand_specific"] == 0 ||
                json["isoform_parameters"]["strand_specific"] == 1) &&
            json["isoform_parameters"]["remove_incomp_reads"] >= 0)) {
        return 0;
    }

    // check that these values are valid booleans
    if (!json["global_parameters"]["generate_raw_isoform"].isBool() ||
        !json["global_parameters"]["has_UMI"].isBool()) {
        return 0;
    }

    // the file is valid
    return 1;
}

Config load_json_config(Json::Value json) {
    /*
        takes a root json value object,
        reads all the values and populates a config object
    */

    Config config;

    config.pipeline_parameters.do_genome_alignment = json["pipeline_parameters"]["do_genome_alignment"].asBool();
    config.pipeline_parameters.do_isoform_identification = json["pipeline_parameters"]["do_isoform_identification"].asBool();
    config.pipeline_parameters.do_read_realignment = json["pipeline_parameters"]["do_read_realignment"].asBool();
    config.pipeline_parameters.do_transcript_quantification = json["pipeline_parameters"]["do_transcript_quantification"].asBool();

    config.global_parameters.generate_raw_isoform = json["global_parameters"]["generate_raw_isoform"].asBool();
    config.global_parameters.has_UMI = json["global_parameters"]["has_UMI"].asBool();

    config.isoform_parameters.MAX_DIST = json["isoform_parameters"]["MAX_DIST"].asInt();
    config.isoform_parameters.MAX_TS_DIST = json["isoform_parameters"]["MAX_TS_DIST"].asInt();
    config.isoform_parameters.MAX_SPLICE_MATCH_DIST = json["isoform_parameters"]["MAX_SPLICE_MATCH_DIST"].asInt();
    config.isoform_parameters.MIN_FL_EXON_LEN = json["isoform_parameters"]["min_fl_exon_len"].asInt();
    config.isoform_parameters.MAX_SITE_PER_SPLICE = json["isoform_parameters"]["Max_site_per_splice"].asInt();
    config.isoform_parameters.MIN_SUP_CNT = json["isoform_parameters"]["Min_sup_cnt"].asInt();

    config.isoform_parameters.MIN_CNT_PCT = json["isoform_parameters"]["Min_cnt_pct"].asFloat();
    config.isoform_parameters.MIN_SUP_PCT = json["isoform_parameters"]["Min_sup_pct"].asFloat();
    config.isoform_parameters.STRAND_SPECIFIC = json["isoform_parameters"]["strand_specific"].asInt();
    config.isoform_parameters.REMOVE_INCOMP_READS = json["isoform_parameters"]["remove_incomp_reads"].asInt();

    config.alignment_parameters.use_junctions = json["alignment_parameters"]["use_junctions"].asBool();
    config.alignment_parameters.no_flank = json["alignment_parameters"]["no_flank"].asBool();

    config.realign_parameters.use_annotation = json["realign_parameters"]["use_annotation"].asBool();

    config.transcript_counting.min_tr_coverage = json["transcript_counting"]["min_tr_coverage"].asFloat();
    config.transcript_counting.min_read_coverage = json["transcript_counting"]["min_read_coverage"].asFloat();

    return config;
}