#ifndef PARSING_H
#define PARSING_H

#include <string>
#include <vector>
#include <cctype>
#include <sstream>
#include <fstream>

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

inline std::string leftStrip(const std::string &str) {
    size_t firstCharPos = 0;
    while (isspace(str[firstCharPos])) firstCharPos++;
    return str.substr(firstCharPos);
}

inline std::string rightStrip(const std::string &str) {
    size_t firstCharPos = str.size() - 1;
    while (isspace(str[firstCharPos])) firstCharPos--;
    return str.substr(0, firstCharPos + 1);
}

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
#endif // PARSING_H