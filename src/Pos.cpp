#include "Pos.h"

// Pos::Pos(std::string chr, int start, int end, char strand, std::string parent_id)
// {
//     this->chr = chr;
//     this->start = start;
//     this->end = end;
//     this->strand = strand;
//     this->parent_id = parent_id;
// }

// Pos::Pos(i)
// {
//     this->chr = "aaaa";
// }

// using namespace Rcpp;

// Rcpp::List
// Pos::to_R()
// {
//     /*
//         wraps up the Pos object into an Rcpp list
//     */
    
//     return List::create(
//         _["chr"] = this->chr,
//         _["start"] = this->start,
//         _["end"] = this->end,
//         _["strand"] = this->strand,
//         _["parent_id"] = this->parent_id
//     );
// }

using namespace Rcpp;

List
pos_to_R(Pos * pos)
{
    /*
        wraps up the Pos struct into an Rcpp list
    */
    
    return List::create(
        _["chr"] = pos->chr,
        _["start"] = pos->start,
        _["end"] = pos->end,
        _["strand"] = pos->strand,
        _["parent_id"] = pos->parent_id
    );
}

Pos
pos_from_R(List list)
{
    Pos pos;

    // pos.chr = list["chr"];
    // pos.start = 
    // pos.parent_id = list["parent_id"];
    return pos;
}

// void
// Pos::from_R(List list)
// {
//     /*
//         unwraps an Rcpp list and imports all the values
//     */

//     this->chr = std::string(list["chr"]);
//     this->start = atoi(list["start"]);
//     this->end = list["end"];
//     this->strand = list["strand"];
//     this->parent_id = std::string(list["parent_id"]);
// }

// void
// Pos::print()
// {
//     /*
//         prints out the data of the Pos object for debugging
//     */

//     std::cout << "Pos\n"
//             << "\tchr : " << this->chr << "\n"
//             << "\tstart : " << this->start << "\n"
//             << "\tend : " << this->end << "\n"
//             << "\tstrand : " << this->strand << "\n"
//             << "\tparent_id : " << this->parent_id << "\n";
// }