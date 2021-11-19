#include <string>
#include <map>
#include <vector>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <Rcpp.h>

#include "GeneBlocks.h"
#include "junctions.h"
#include "misc.h"
#include "Isoforms.h"
#include "bam.h"
#include "StartEndPair.hpp"
#include "cigars.h"
#include "BamRecord.hpp"

using namespace Rcpp;

#ifndef ISOKEY
#define ISOKEY

/*  struct that we can use as a key in a map
*/
struct IsoformKey {
    std::string chr;
    int start;
    int end;

    // we need this to correctly logically compare IsoformKeys
    bool operator==(const IsoformKey &other) const
    { 
        return ((chr == other.chr) &&
                (start == other.start) && 
                (end == other.end));
    }

    bool operator<(const IsoformKey &other) const
    {
        // compare a and b, return true if a is 'less than' b
        if (start < other.start) {
            return true;
        } else if ((start == other.start) && (end < other.end)) {
            return true;
        }
        return false;
    }
};


namespace std {
    template <> struct hash<IsoformKey>
    {
        std::size_t operator()(const IsoformKey& k) const
        {
            using std::size_t;
            using std::hash;

            return (hash<string>()(k.chr))
                ^ ((hash<int>()(k.start)
                ^ (hash<int>()(k.end) << 1)) >> 1);
        }
    };
}

/*  quick struct so we have something to pass down both ref name and records
    to the BAM fetch_function
*/
struct DataStruct {
    bam_header_t * header;
    std::vector<BAMRecord> * records;
};

#define BAM_CMATCH      0   // CIGAR character for matching
#define BAM_CDEL        2   
#define BAM_CREF_SKIP   3
#define BAM_CEQUAL      7
#define BAM_CDIFF       8

#endif

static int
fetch_function(const bam1_t *b, void *data);

void
bam_read (std::string bam_in, std::string chr, int s, int e);

std::vector<StartEndPair>
get_blocks(BAMRecord record);

void
group_bam2isoform (
    std::string bam_in, 
    std::string out_gff3, 
    std::string out_stat, 
    std::unordered_map<std::string, std::vector<GeneBlocks>>    * chr_to_blocks, 
    std::unordered_map<std::string, std::vector<StartEndPair>>  * gene_dict, 
    std::unordered_map<std::string, Junctions>                  * transcript_to_junctions,
    std::unordered_map<std::string, Pos>                        * transcript_dict,
    std::string fa_f,
    IsoformParameters isoform_parameters,
    std::string raw_gff3
);