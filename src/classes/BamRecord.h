#ifndef BAM_RECORD_H
#define BAM_RECORD_H

#include <vector>
#include <string>

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

	inline bool operator==(const Flag &o) const{
		return
			read_paired == o.read_paired &&
			read_mapped_in_proper_pair == o.read_mapped_in_proper_pair &&
			read_unmapped == o.read_unmapped &&
			mate_unmapped == o.mate_unmapped &&
			read_reverse_strand == o.read_reverse_strand &&
			mate_reverse_strand == o.mate_reverse_strand  &&
			first_in_pair == o.first_in_pair &&
			second_in_pair == o.second_in_pair &&
			not_primary_alignment == o.not_primary_alignment &&
			read_fails_platform_vendor_quality_checks == o.read_fails_platform_vendor_quality_checks &&
			read_is_PCR_or_optical_duplicate == o.read_is_PCR_or_optical_duplicate &&
			supplementary_alignment == o.supplementary_alignment;
	}
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


Flag
read_flag(int);

BAMRecord
read_record(const bam1_t*, const bam_header_t*);

int calculate_length_from_cigar(const std::vector<CigarPair> &cigartuples, bool queryLength);

#endif // BAM_RECORD_H
