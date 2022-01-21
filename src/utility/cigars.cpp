/*
  All the code relating to cigars is in here
*/

#include "cigars.h"

#include <string>
#include <vector>

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
generate_cigar (std::vector<CigarPair> cigar)
{
  /*
    takes a cigar as a vector of pairs of ints,
    converts it to CIGAR string format
  */

  char CIGAR_CODE[] = "MIDNSHP=XB";
  std::string result = "";

  for (const auto & pair : cigar) {
    result.append(std::to_string(pair.len));
    result.push_back(CIGAR_CODE[pair.op]);
  }
  return result;
}

std::vector<CigarPair>
smooth_cigar (std::vector<CigarPair> cigar, int threshold)
{
  /*
    takes a cigar as a vector of int pairs,
    smooths it down,
    returns the smooth cigar as a vector of int pairs
  */

  std::vector<CigarPair> new_cigar = {cigar[0]};
  
  for (int i = 1; i < cigar.size(); i++) {
    if (new_cigar.back().op != 0) {
      new_cigar.push_back(cigar[i]);
    } else if (cigar[i].op == 0) { // merge matched reads 
      new_cigar.back().len = new_cigar.back().len + cigar[i].len;
    } else if (cigar[i].op == 2) {
      if (cigar[i].len <= threshold) {
        new_cigar.back().len = new_cigar.back().len + cigar[i].len;
      } else {
        new_cigar.push_back(cigar[i]);
      }
    } else if (cigar[i].op == 1) {
      if (cigar[i].len > threshold) {
        new_cigar.push_back(cigar[i]);
      } else {
        continue;
      }
    } else {
      new_cigar.push_back(cigar[i]);
    }
  }
  return new_cigar;
}
