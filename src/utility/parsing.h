#ifndef PARSING_H
#define PARSING_H

#include <string>
#include <vector>
#include <cctype>
#include <sstream>

inline std::vector<std::string> splitStringToVector(const std::string &str, const char sep) {
    std::vector<std::string> tokens;
    std::stringstream ss(str);
    std::string line;
    while (std::getline(ss, line, sep)) {
        tokens.push_back(line);
    }
    return tokens;
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
#endif // PARSING_H