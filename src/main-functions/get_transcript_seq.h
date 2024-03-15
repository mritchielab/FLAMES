#ifndef GET_TRANSCRIPT_SEQ_H
#define GET_TRANSCRIPT_SEQ_H

#include <string>

#include "../classes/GFFData.h"

void
get_transcript_seq(
    const std::string &fa_file,
    const std::string &fa_out,
    const GFFData &isoform_annotation,
    const GFFData &ref_annotation
);

#endif // GET_TRANSCRIPT_SEQ_H