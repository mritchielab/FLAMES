#ifndef PARSING_H
#define PARSING_H

#include <Rcpp.h>
#include <string>
#include <vector>
#include <cctype>
#include <sstream>
#include <fstream>
#include <stdexcept> // invalid_argument

inline std::vector<std::string> splitStringToVector(const std::string &str, const char sep) {
    std::vector<std::string> tokens;
    std::stringstream ss(str);
    std::string line;
    while (std::getline(ss, line, sep)) {
        tokens.push_back(line);
    }
    return tokens;
}
inline std::string getFilenameBeforeExt(const std::string &str, const char sep='.') {
    int i = str.size() - 1;
    while (str[i] != sep) i--;
    return str.substr(0, i);
}

// Remove whitespace characters from the left side of a string
inline std::string leftStrip(const std::string &str) {
    size_t firstCharPos = 0;
    while (isspace(str[firstCharPos])) firstCharPos++;
    return str.substr(firstCharPos);
}

// Remove whitespace characters from the right of a string
inline std::string rightStrip(const std::string &str) {
    size_t firstCharPos = str.size() - 1;
    while (isspace(str[firstCharPos])) firstCharPos--;
    return str.substr(0, firstCharPos + 1);
}

// Remove whitespace characters from both sides of a string
// return: a new string without whitespace on left and right sides
inline std::string strip(const std::string &str) {
    std::string intermediate = leftStrip(str);
    return rightStrip(intermediate);
}

inline void foreachLineinFile(const std::string &filename, const std::function<void(const std::string &)> func) {
    std::ifstream infile(filename.c_str());
    if (!infile.is_open()) return;

    std::string line;
    while (std::getline(infile, line)) func(line);
    infile.close();
}

inline int parseInt(const std::string &s, size_t *idx = 0, int base = 10) {
    try {
        return std::stoi(s, idx, base);
    } catch (const std::invalid_argument  &err) {
        Rcpp::Rcout << "invalid argument to stoi: " << s << ". Using -1 instead.\n";
        return -1;
    }
    return 0;
}
inline float parseFloat(const std::string &s, size_t *idx = 0) {
    try {
        return std::stof(s, idx);
    } catch (const std::invalid_argument  &err) {
        Rcpp::Rcout << "invalid argument to stoi: " << s << ". Using -1 instead.\n";
        return -1.0;
    }
    return 0.0;
}
#endif // PARSING_H