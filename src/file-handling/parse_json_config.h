#ifndef PARSE_JSON_CONFIG_H
#define PARSE_JSON_CONFIG_H

#include <string>
#include <Rcpp.h>

#include "json/json.h"
#include "config.h"

Rcpp::List
parse_json_config_cpp (std::string);

int
verify_json_config(Json::Value);

Config
load_json_config(Json::Value);

void print_config_cpp(Rcpp::List);

#endif // PARSE_JSON_CONFIG_H