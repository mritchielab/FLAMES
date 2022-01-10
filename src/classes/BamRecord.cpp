#include "BamRecord.h"

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
    const uint32_t *cigar = bam_get_cigar(b);
    for (int k = 0; k < b->core.n_cigar; k++) {
        cigar_pairs.push_back((CigarPair){
            (int)bam_cigar_op(cigar[k]),
            (int)bam_cigar_oplen(cigar[k])
        });
    }
    return cigar_pairs;
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
    rec.reference_start = b->core.pos;
    rec.reference_end = b->core.pos + b->core.l_qseq;
    rec.reference_name = header->target_name[b->core.tid];
    rec.AS_tag = bam_aux2i(bam_aux_get(b, "AS"));
    rec.query_alignment_length = b->core.l_qseq;
    rec.read_name = bam_get_qname(b);
    rec.mapping_quality = (int)(b->core.qual);

    // rec.reference_name = std::string(sam_hdr_tid2name(header, b->core.tid));
    // std::cout << sam_hdr _tid2len(header, b->core.tid) << "\n";
    rec.cigar = generate_cigar_pairs(b);
    rec.flag = read_flag(b->core.flag);
    return rec;
}