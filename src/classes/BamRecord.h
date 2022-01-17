#ifndef RECORD
#define RECORD

#include <vector>
#include <string>
#include <codecvt>
#include <Rcpp.h>

#include "../utility/bam.h"
#include "../utility/cigars.h"

/* 
    holds all of the flags in a bam record
*/
struct Flag {
    bool
    read_paired;
    bool
    read_mapped_in_proper_pair;
    bool
    read_unmapped;
    bool
    mate_unmapped;
    bool
    read_reverse_strand;
    bool
    mate_reverse_strand;
    bool
    first_in_pair;
    bool
    second_in_pair;
    bool
    not_primary_alignment;
    bool
    read_fails_platform_vendor_quality_checks;
    bool
    read_is_PCR_or_optical_duplicate;
    bool
    supplementary_alignment;
};

/* 
    a struct for handling each record in a bam file
*/
struct BAMRecord {
    std::vector<CigarPair>
    cigar;
    std::string
    cigar_string;

    std::string
    reference_name;
    int
    reference_start;
    int
    reference_end;
    
    std::string
    read_name;

    float
    AS_tag;
    int
    query_alignment_length;
    int
    mapping_quality;
    Flag
    flag;
};

#endif

std::vector<CigarPair>
generate_cigar_pairs(const bam1_t*);

Flag
read_flag(int);

BAMRecord
read_record(const bam1_t*, const bam_header_t*);