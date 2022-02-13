#ifndef ISOKEY
#define ISOKEY

#include <string>

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

#endif // ISOKEY