#ifndef FASTQ_UTILS_H
#define FASTQ_UTILS_H

#include <string>

#include <Rcpp.h>
#include <R.h>
#include "zlib.h"
#include "htslib/kseq.h"

#ifndef INIT_KSEQ
#define INIT_KSEQ
KSEQ_INIT(gzFile, gzread)
#endif

void REMOVE_WARNINGS();

void fq_gz_write(gzFile out_file, std::string name, std::string qual, std::string seq);

void fq_gz_write(gzFile out_file, kseq_t *seq);
#endif // #ifndef FASTQ_UTILS_H