#include <iostream>
#include <string>
#include <vector>


#ifndef CIG_PAIR
#define CIG_PAIR

/* a struct for handling each CIGAR operation
*/
struct CigarPair
{
    int op;     // which operation to perform, ranged 0 to 4
    int len;    // length of operation
};

#endif

std::string
generate_cigar (std::vector<CigarPair> cigar);

std::vector<CigarPair>
smooth_cigar (std::vector<CigarPair> cigar, int threshold);