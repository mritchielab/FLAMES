#include <Rcpp.h>
#include <R.h>
#include "zlib.h"
#include "kseq.h"
#include <iostream>
#include <sstream>
#include <fstream>

using namespace Rcpp;

#ifndef INIT_KSEQ
#define INIT_KSEQ
KSEQ_INIT(gzFile, gzread)
#endif

#ifndef FQ_GZ_WRITE
#define FQ_GZ_WRITE
void fq_gz_write(gzFile out_file, std::string name, std::string qual, std::string seq);

void fq_gz_write(gzFile out_file, kseq_t *seq);
#endif