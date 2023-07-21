#ifndef FIND_ISOFORM_H
#define FIND_ISOFORM_H

#include <Rcpp.h>

#include <string>

#include "../classes/Pos.h"

void
find_isoform_multithread_cpp
(
    const std::string &gff3,
    const std::string &genome_bam,
    const std::string &isoform_gff3,
    const std::string &tss_tes_stat,
    const std::string &genomefa,
    const std::string &transcript_fa,
    const Rcpp::List  &isoform_parameters,
    const std::string &raw_splice_isoform);

#endif // FIND_ISOFORM_H