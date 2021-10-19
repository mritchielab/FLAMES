#ifndef POS
#define POS

#include <string>
#include <Rcpp.h>
#include <iostream>

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

Rcpp::List
pos_to_R(Pos * pos);

Pos
pos_from_R(Rcpp::List list);

#endif