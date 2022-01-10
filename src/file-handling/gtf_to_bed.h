#ifndef GTF_TO_BED_H
#define GTF_TO_BED_H

#include <string>

void
gtf_to_bed_cpp(std::string in_gtf, std::string out_bed, std::string chrom_sizes_file);

#endif // GTF_TO_BED_H