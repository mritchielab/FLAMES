#ifndef TYPES_H
#define TYPES_H
/*
Header file utility types used to improve readability of
code.
exon: An exon is a vector of StartEndPairs
transcriptvector: a vector of transcript strings
*/

#include <string>
#include <vector>

#include "StartEndPair.h"

typedef StartEndPair exon;

typedef std::vector<std::string> transcriptvector;

#endif // TYPES_H