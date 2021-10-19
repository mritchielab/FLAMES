#ifndef START_END_PAIR
#define START_END_PAIR

#include <vector>

struct StartEndPair {
    int start;
    int end;

    // we need this to correctly logically compare StartEndPairs
    bool operator==(const StartEndPair &other) const
    { 
        return (start == other.start && end == other.end);
    }

    bool operator<(const StartEndPair &other) const
    {
        // compare a and b, return true if a is 'less than' b
        // in this case, 'less than' is defined if a.start is less than b.start
        return (start < other.start);
    }
};

/*  next, we will need hashing functions
    this is so that we can use a std::vector<StartEndPair> as a key in a dictionary
*/
namespace std {
    template <>
        struct hash<StartEndPair>
        {
            std::size_t operator()(const StartEndPair& k) const
            {
                using std::size_t;
                using std::hash;

                return ((hash<int>()(k.start)
                    ^ (hash<int>()(k.end) << 1)) >> 1);
            }
        };

    template<>
        struct hash<vector<StartEndPair>>
        {
            std::size_t operator()(const vector<StartEndPair>& vec) const
            {
                using std::size_t;
                using std::hash;

                std::size_t seed = vec.size();

                for (auto& pair : vec) {
                    seed ^= ((hash<int>()(pair.start) ^ (hash<int>()(pair.end))) >> 1);
                }
                return seed;
            }
        };
}

#endif