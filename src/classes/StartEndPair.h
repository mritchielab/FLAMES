#ifndef START_END_PAIR_H
#define START_END_PAIR_H

#include <vector>
#include <sstream>
#include <string>

struct StartEndPair {
    int start;
    int end;

    StartEndPair(int _start, int _end): start(_start), end(_end) {}
    StartEndPair(const StartEndPair &sep): start{sep.start}, end{sep.end} {}

    // we need this to correctly logically compare StartEndPairs
    bool operator==(const StartEndPair &other) const
    { 
        return (start == other.start && end == other.end);
    }
	bool operator!=(const StartEndPair &other) const
	{
		return (start != other.start || end != other.end);
	}

    bool operator<(const StartEndPair &other) const
    {
        // compare a and b, return true if a is 'less than' b
        if (start < other.start) {
            return true;
        } else if ((start == other.start) && (end < other.end)) {
            return true;
        }
        return false;
    }

    bool operator>(const StartEndPair &other) const
    {
        // compare a and b, return true if a is 'greater than' b
        if (start > other.start) {
            return true;
        } else if ((start == other.start) && (end > other.end)) {
            return true;
        }
        return false;
    }

    bool operator>=(const StartEndPair &other) const
    {
        return ((*this) > other) || ((*this) == other);
    }

    bool operator<=(const StartEndPair &other) const
    {
        return ((*this) < other) || ((*this) == other);
    }

    std::string getString() const {
        std::stringstream s;
        s << "(" << start << ", " << end << ")";
        return s.str();
    }
};

/*
    next, we will need hashing functions
    this is so that we can use a std::vector<StartEndPair> as a key in a dictionary
*/
namespace std {
    template <> struct hash<StartEndPair>
    {
        std::size_t operator()(const StartEndPair& k) const
        {
            using std::size_t;
            using std::hash;

            return ((hash<int>()(k.start)
                ^ (hash<int>()(k.end) << 1)) >> 1);
        }
    };

    template<> struct hash<vector<StartEndPair>>
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
    
    template <> struct hash<vector<int>>
    {
        size_t operator()(vector<int> const& vec) const 
        {
            size_t seed = vec.size();
            for(auto& i : vec) {
                seed ^= i + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            }
            return seed;
        }
    };
}

inline bool StartEndPairCompare(const StartEndPair &a, const StartEndPair &b) {
    // compare a and b, return true if a is 'less than' b
    // in this case, 'less than' is defined if a.start is less than b.start
    return a.start < b.start;
}

#endif // START_END_PAIR_H