#include <string>
#include <Rcpp.h>

#ifndef POS
#define POS

class Pos
{
    public:
        std::string chr;
        int start;
        int end;
        char strand;
        std::string parent_id;

        Rcpp::List
        to_R();
        void
        from_R(Rcpp::List list);
        void
        print();
};

#endif