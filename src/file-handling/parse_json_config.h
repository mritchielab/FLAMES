#ifndef PARSE_JSON_CONFIG_H
#define PARSE_JSON_CONFIG_H

#include <string>
#include <Rcpp.h>

#include "../utility/json/json.h"
#include "../classes/Config.h"

Rcpp::List
parse_json_config(std::string);

int
verify_json_config(Json::Value);

Config
load_json_config(Json::Value);

void print_config(Rcpp::List);

#endif // PARSE_JSON_CONFIG_H