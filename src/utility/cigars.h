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

typedef std::vector<CigarPair> Cigar; 

#endif

std::string
generate_cigar (Cigar cigar);

Cigar
smooth_cigar (Cigar cigar, int threshold);