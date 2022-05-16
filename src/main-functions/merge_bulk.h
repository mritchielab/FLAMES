#ifndef MERGE_BULK_H
#define MERGE_BULK_H

#include <Rcpp.h>
#include <R.h>
#include "zlib.h"
#include "htslib/kseq.h"
#include <iostream>
#include <sstream>
#include <fstream>

#include "../utility/fastq_utils.h"

const char * 
shorten_filename
(
    const char *file_name, 
    int length, 
    int &out_length
);

void
merge_bulk_fastq
(
    Rcpp::StringVector fastq_files,
    Rcpp::String out_fastq
);

#endif