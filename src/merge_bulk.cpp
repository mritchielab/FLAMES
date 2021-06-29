#include "fastq_utils.h"

using namespace Rcpp;


//' Merge Bulk Fastq Files
//' 
//' @description Merge all fastq files into a single fastq
//' This function also inserts the original fastq file name into the header line of each read.
//'
//' @param fastq_files a string vector of the fastq file paths to merge
//' @param out_fastq the fastq file path of the output file
//' @return returns NULL
//' @useDynLib FLAMES, .registration=TRUE
// [[Rcpp::export]]
void merge_bulk_fastq_cpp(StringVector fastq_files, String out_fastq) {
    gzFile fp;
    kseq_t *seq;
    int l;

    gzFile o_stream_gz = gzopen(out_fastq.get_cstring(), "wb2");
    
    const char *separator = "_NNN#";
    for (int i = 0; i < fastq_files.size(); i++) {
        // For every fastq file, read in each read and prefix the name line with the file name
        fp = gzopen(fastq_files(i), "r");
        seq = kseq_init(fp);

        String file_name = fastq_files(i);
        const char *c_file_name = file_name.get_cstring();
        int file_name_length = fastq_files(i).size();
        int offset = file_name_length + 5;
        
        while ((l = kseq_read(seq)) >= 0) {
            // reallocate name block to expand for {filename}_NNN#{seq->name.s}
            seq->name.s = (char *)realloc(seq->name.s, offset + seq->name.l);

            // move name along, and insert file name
            char *const seq_name = seq->name.s;
            memmove(seq_name + offset, seq_name, (seq->name.l + 1) * sizeof(char));
            memcpy(seq_name, c_file_name, file_name_length);
            memcpy(seq_name + file_name_length, separator, 5);

            fq_gz_write(o_stream_gz, seq);
        }

        kseq_destroy(seq);
        gzclose(fp);
    }

    gzclose(o_stream_gz);
}