#ifndef POS_H
#define POS_H

#include <string>
#include <Rcpp.h>

struct Pos 
{
    std::string chr;
    int start;
    int end;
    char strand;
    std::string parent_id;

    Pos() {}
    Pos(std::string _chr, int _start, int _end, char _strand, std::string _parent_id)
        : chr(_chr), start(_start), end(_end), strand(_strand), parent_id(_parent_id) {}
};

inline bool comparePos(const Pos &a, const Pos &b) {
	return a.chr == b.chr 
		&& a.start == b.start
		&& a.end == b.end
		&& a.strand == b.strand
		&& a.parent_id == b.parent_id;
}

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
    Pos pos(
        (Rcpp::String)(list["chr"]),
        list["start"],
        list["end"],
        list["strand"],
        (Rcpp::String)list["parent_id"]
    );
    return pos;
}

#endif // POS_H