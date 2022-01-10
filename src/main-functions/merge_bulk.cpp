#include <Rcpp.h>
#include <R.h>
#include "zlib.h"
#include "htslib/kseq.h"
#include <iostream>
#include <sstream>
#include <fstream>

#include "../utility/fastq_utils.h"

using namespace Rcpp;

const char * shorten_filename(const char *file_name, int length, int &out_length) {
    int slash = -1;
    for (int i = length - 1; i >= 0; i--) { 
        if (file_name[i] == '/') {
            slash = i;
            break;
        }   
    }
    const char *short_name = file_name + (slash + 1);
    out_length = length - slash - 1;
    return short_name;
}

// This has the potential to be parralellized if time allows.
// create a thread per file, and process to the one out file using a file handler object with mutex lock on file
//' Merge Bulk Fastq Files
//' 
//' @description Merge all fastq files into a single fastq
//' This function also inserts the original fastq file name into the header line of each read.
//'
//' @param fastq_files a string vector of the fastq file paths to merge
//' @param out_fastq the fastq file path of the output file
//' @return returns NULL
//' @useDynLib FLAMES, .registration=TRUE
//' @import zlibbioc
// [[Rcpp::export]]
void merge_bulk_fastq_cpp(StringVector fastq_files, String out_fastq) {
    gzFile fp;
    kseq_t *seq;
    int l;

    gzFile o_stream_gz = gzopen(out_fastq.get_cstring(), "wb2");
    Rcout << o_stream_gz << "\n";
    
    // int array to track the number of reads processed in each fastq file.
    unsigned int * read_counts = (unsigned int *)malloc(fastq_files.size() * sizeof(unsigned));

    const char *separator = "_NNN#"; // separator string between new read name and old read name

    for (short unsigned int i = 0; i < fastq_files.size(); i++) {
        // For every fastq file, read in each read and prefix the name line with the file name
        fp = gzopen(fastq_files(i), "r");
        seq = kseq_init(fp);

        read_counts[i] = 0;

        String file_name = fastq_files(i);
        const char *c_file_name = file_name.get_cstring();
        int file_name_length = fastq_files(i).size();
        //shorten the file name to only include local name (not full path name)
        c_file_name = shorten_filename(c_file_name, file_name_length, file_name_length);
        // c_file_name and file_name_length are now both for the shortened versions

        int offset = file_name_length + 5;
        
        while ((l = kseq_read(seq)) >= 0) {
            read_counts[i]++;
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

        Rcout << c_file_name << ": " << read_counts[i] << "\n";
    }

    free(read_counts);
    gzclose(o_stream_gz);
}