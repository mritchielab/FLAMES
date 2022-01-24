#ifndef ISOKEY
#define ISOKEY

#include <string>
#include <vector>

#include "../classes/BamRecord.h"
#include "../utility/bam.h"

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

#endif // ISOKEY