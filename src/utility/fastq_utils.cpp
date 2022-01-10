#include "fastq_utils.hpp"


#ifndef INIT_KSEQ
#define INIT_KSEQ
KSEQ_INIT(gzFile, gzread)
#endif

void REMOVE_WARNINGS() {
    (void)&kseq_init;
    (void)&kseq_read;
    (void)&kseq_destroy;
}

void fq_gz_write(gzFile out_file, std::string name, std::string qual, std::string seq) {
    (void)REMOVE_WARNINGS();
    std::stringstream stream;
    stream << "@" << name << "\n" <<
        seq << "\n" <<
        "+" << "\n" <<
        qual << "\n";
    gzputs(out_file, stream.str().c_str());
}

void fq_gz_write(gzFile out_file, kseq_t *seq) {
    std::stringstream stream;
    stream << "@" << seq->name.s << "\n" <<
        (seq->seq.s) << "\n" <<
        "+" << "\n" <<
        (seq->qual.s) << "\n";
    gzputs(out_file, stream.str().c_str());
}