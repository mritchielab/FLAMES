#include "cigars.h"

#include <string>
#include <vector>
#include <sstream>

#include "bam.h"

/*  take a bam entry,
    populate a vector with all of its CIGAR operations
    this is to mimic the output of bamnostic's default cigar from bamfile.fetch
*/
std::vector<CigarPair>
generate_cigar_pairs(const bam1_t *b)
{
    std::vector<CigarPair>
    cigar_pairs;

    // iterate over the cigar
    const uint32_t *cigar = bam1_cigar(b);
    for (int k = 0; k < (int)bam1_cigar_len(b); k++) {
        cigar_pairs.push_back((CigarPair){
            (int)bam_cigar_op(cigar[k]),
            (int)bam_cigar_oplen(cigar[k])
        });
    }
    return cigar_pairs;
}

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

std::string
generate_cigar (const std::vector<CigarPair> &cigar)
{
  /*
    takes a cigar as a vector of pairs of ints,
    converts it to CIGAR string format
  */

  const char CIGAR_CODE[] = "MIDNSHP=XB";
  std::string result = "";

  for (const auto & pair : cigar) {
    result.append(std::to_string(pair.len));
    result.push_back(CIGAR_CODE[pair.op]);
  }
  return result;
}

/*
takes a cigar as a vector of int pairs,
smooths it down,
returns the smooth cigar as a vector of int pairs
*/
std::vector<CigarPair>
smooth_cigar (const std::vector<CigarPair> &cigar, int threshold)
{
    std::vector<CigarPair> new_cigar = {cigar[0]};

    for (int i = 1; i < (int)cigar.size(); i++) {
        if (new_cigar.back().op != 0) {
            new_cigar.push_back(cigar[i]);
            continue;
        } 
        
        switch(cigar[i].op) {
            case 0:
                new_cigar.back().len += cigar[i].len;
                break;
            case 1:
                if (cigar[i].len > threshold) {
                    new_cigar.push_back(cigar[i]);
                } 
                break;
            case 2:
                if (cigar[i].len <= threshold) {
                    new_cigar.back().len += cigar[i].len;
                } else {
                    new_cigar.push_back(cigar[i]);
                }
                break;
            default:
                new_cigar.push_back(cigar[i]);
        }
    }

    return new_cigar;
}

std::string printCigarPairs(const std::vector<CigarPair> &cigar) {
    std::stringstream ss;
    ss << cigar[0].len << "-" << cigar[0].op;
    for (int i =1; i < cigar.size(); i++) {
        ss << "," << cigar[i].len << "-" << cigar[i].op;
    }
    return ss.str();
}