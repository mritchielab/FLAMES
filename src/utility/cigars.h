#ifndef CIG_PAIR
#define CIG_PAIR

#include <string>
#include <vector>

#include "bam.h"

/*
cigar specification:
(lifted from the original sc_longread.py)

M   BAM_CMATCH  0
I   BAM_CINS    1
D   BAM_CDEL    2
N   BAM_CREF_SKIP   3
S   BAM_CSOFT_CLIP  4
H   BAM_CHARD_CLIP  5
P   BAM_CPAD    6
=   BAM_CEQUAL  7
X   BAM_CDIFF   8
B   BAM_CBACK   9

*/

/* a struct for handling each CIGAR operation
*/
struct CigarPair
{
    int op;     // which operation to perform, ranged 0 to 4
    int len;    // length of operation

	inline bool operator==(const CigarPair &a) const {
		return (a.op == op) && (a.len == len);
	}
};

std::vector<CigarPair>
generate_cigar_pairs(const bam1_t*);

std::string
generate_cigar (std::vector<CigarPair> cigar);

std::vector<CigarPair>
smooth_cigar (std::vector<CigarPair> cigar, int threshold);

#endif