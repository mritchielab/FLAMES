#ifndef DATA_STRUCT_H
#define DATA_STRUCT_H

#include <vector>

#include "../utility/bam.h"
#include "BamRecord.h"

/*  quick struct so we have something to pass down both ref name and records
    to the BAM fetch_function
*/
struct DataStruct {
    bam_header_t * header;
    std::vector<BAMRecord> * records;
};

#define BAM_CMATCH      0   // CIGAR character for matching
#define BAM_CDEL        2   
#define BAM_CREF_SKIP   3
#define BAM_CEQUAL      7
#define BAM_CDIFF       8

#endif // DATA_STRUCT_H
