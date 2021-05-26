#ifndef VECTORCOUNTER_H
#define VECTORCOUNTER_H

#include <Rcpp.h>
using namespace Rcpp;

class VectorEnumerator {
        std::vector<int> orders;
        std::vector<NumericVector> cases;
    public:
        VectorEnumerator(int maxOrder);
        VectorEnumerator& extendDimension(int coordinateMaxOrder, int maxOrder);
        NumericMatrix getAllCases();
};

class VectorCounter {
        int count;
        int maxcount;
        NumericMatrix table;
    public:
        VectorCounter(NumericVector max, int order);
        VectorCounter(NumericVector max);
        int getCount();
        int getMaxcount();
        NumericMatrix getTable();
        NumericVector degree();
        void increase();
        bool inRange();
};

#endif