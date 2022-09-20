#include <vector>
#include <string>
#include <codecvt>
#include <Rcpp.h>

#include "../utility/bam.h"
#include "../utility/cigars.h"
#include "BamRecord.h"

int calculate_length_from_cigar(const std::vector<CigarPair> &cigartuples, bool queryLength) {
	// calculate the length of either query or reference using a cigar string
	// query consuming operations are: MIS=X (0, 1, 4, 7, 8)
	// reference consuming operations are: MDN=X (0, 2, 3, 7, 8)
	// if queryLength =True, calculate query length, otherwise calculate read length
	int len = 0;
	for (const auto &pair : cigartuples) {
		if (queryLength) {
			if (pair.op == 0 || pair.op == 1 || pair.op == 4 || pair.op == 7 || pair.op == 8) {
				len += pair.len;
			}
		} else {
			if (pair.op == 0 || pair.op == 2 || pair.op == 3 || pair.op == 7 || pair.op == 8) {
				len += pair.len;
			}
		}
	}

	return len;
}

int calculate_reference_end(const bam1_t *b, const std::vector<CigarPair> &cigartuples) {
	// calculate the length of the reference based on reference consuming operations:
	// MDN=X (0, 2, 3, 7, 8)
	int ref_start = bam_reference_start(b);
	ref_start += calculate_length_from_cigar(cigartuples, false);
	return ref_start;
}

int calculate_query_alignment_length(const std::vector<CigarPair> &cigartuples) {
	int q_len = calculate_length_from_cigar(cigartuples, true);

	// remove regions of soft clipping from the ends of the query length 
	if (cigartuples.size() > 0) {
		if (cigartuples[0].op == 4) {
			q_len -= cigartuples[0].len;
		}
		if (cigartuples.size() > 1) {
			if (cigartuples.back().op == 4) {
				q_len -= cigartuples.back().len;
			}
		}
	}

	return q_len;
}

/*  takes a flag int, converts it to all of the properties it encodes for 
*/
Flag
read_flag(int n)
{
    Flag flag;

    flag.read_paired                               = (n & (1<<0));
    flag.read_mapped_in_proper_pair                = (n & (1<<1));
    flag.read_unmapped                             = (n & (1<<2));
    flag.mate_unmapped                             = (n & (1<<3));
    flag.read_reverse_strand                       = (n & (1<<4));
    flag.mate_reverse_strand                       = (n & (1<<5));
    flag.first_in_pair                             = (n & (1<<6));
    flag.second_in_pair                            = (n & (1<<7));
    flag.not_primary_alignment                     = (n & (1<<8));
    flag.read_fails_platform_vendor_quality_checks = (n & (1<<9));
    flag.read_is_PCR_or_optical_duplicate          = (n & (1<<10));
    flag.supplementary_alignment                   = (n & (1<<11));

    return flag;
}

/*  takes a single bam entry, converts it into a BAMRecord struct
*/
BAMRecord
read_record(const bam1_t * b, const bam_header_t * header)
{
    BAMRecord rec;
    rec.reference_start = bam_reference_start(b);
    rec.reference_name = std::string(header->target_name[b->core.tid]);
    rec.AS_tag = bam_aux2i(bam_aux_get(b, "AS"));
    rec.read_name = std::string(bam1_qname(b));
    rec.mapping_quality = (int)bam_mapping_qual(b);
    rec.cigar = generate_cigar_pairs(b);
	rec.cigar_string = generate_cigar(rec.cigar);
	rec.reference_end = calculate_reference_end(b, rec.cigar);
	rec.query_alignment_length = calculate_query_alignment_length(rec.cigar);
    rec.flag = read_flag(b->core.flag);
    return rec;
}
