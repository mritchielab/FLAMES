#include "Pos.h"
using namespace Rcpp;


List
Pos::to_R()
{
    /*
        wraps up the Pos object into an Rcpp list
    */
    
    return List::create(
        _["chr"] = this->chr,
        _["start"] = this->start,
        _["end"] = this->end,
        _["strand"] = this->strand,
        _["parent_id"] = this->parent_id
    );
}

void
Pos::from_R(List list)
{
    /*
        unwraps an Rcpp list and imports all the values
    */

    this->chr = std::string(list["chr"]);
    this->start = atoi(list["start"]);
    this->end = list["end"];
    this->strand = list["strand"];
    this->parent_id = std::string(list["parent_id"]);
}

void
Pos::print()
{
    /*
        prints out the data of the Pos object
        for debugging mostly
    */

    std::cout << "Pos\n"
            << "\tchr : " << this->chr << "\n"
            << "\tstart : " << this->start << "\n"
            << "\tend : " << this->end << "\n"
            << "\tstrand : " << this->strand << "\n"
            << "\tparent_id : " << this->parent_id << "\n";
}