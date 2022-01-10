#ifndef POS
#define POS

#include <string>
#include <Rcpp.h>

// class Pos 
// {
//     public:
//         Pos (std::string chr, int start, int end, char strand, std::string parent_id);

//         std::string chr;
//         int start;
//         int end;
//         char strand;
//         std::string parent_id;

//         Rcpp::List
//         to_R();
//         void
//         from_R(Rcpp::List list);
//         void
//         print();
// };

struct Pos 
{
    std::string chr;
    int start;
    int end;
    char strand;
    std::string parent_id;
};

inline Rcpp::List pos_to_R(Pos * pos) {
	/*
        wraps up the Pos struct into an Rcpp list
    */
	return Rcpp::List::create(
        Rcpp::_["chr"] = pos->chr,
        Rcpp::_["start"] = pos->start,
        Rcpp::_["end"] = pos->end,
        Rcpp::_["strand"] = pos->strand,
        Rcpp::_["parent_id"] = pos->parent_id
    );
}

inline Pos pos_from_R(Rcpp::List list)
{
    Pos pos;
    pos.chr = (Rcpp::String)(list["chr"]);
    pos.start = list["start"];
    pos.end = list["end"];
    pos.strand = list["strand"];
    pos.parent_id = (Rcpp::String)list["parent_id"];
    return pos;
}

// // [[Rcpp::export]]
// inline void pos_from_R_test(Rcpp::List list)
// {
//     Pos pos = pos_from_R(list);
//     std::cout << "created a pos\n"
//         << "\tchr:" << pos.chr << "\n"
//         << "\tstart:" << pos.start << "\n"
//         << "\tend:" << pos.end << "\n"
//         << "\tstrand:" << pos.strand << "\n"
//         << "\tparent_id:" << pos.parent_id << "\n";
// }

// // [[Rcpp::export]]
// inline Rcpp::List pos_to_R_test()
// {
//     Pos pos = {"chr21", 124, 128, '-', "dan"};
//     std::cout << "created pos on " << pos.chr << "\n";
//     return pos_to_R(&pos);
// }
#endif // POS