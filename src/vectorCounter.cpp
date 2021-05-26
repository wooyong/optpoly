#include "vectorCounter.h"
#include <Rcpp.h>
using namespace Rcpp;

VectorEnumerator::VectorEnumerator(int maxOrder)
{
    // create cases of 0, 1, ..., maxOrder
    for (int k = 0; k <= maxOrder; k++)
    {
        // fill in the members
        orders.push_back(k);
        cases.push_back(NumericVector::create(k));
    }
}

VectorEnumerator &VectorEnumerator::extendDimension(int coordinateMaxOrder, int maxOrder)
{
    // initialize new members
    std::vector<int> newOrders;
    std::vector<NumericVector> newCases;
    // auxiliary variables
    unsigned int i, j, remainingOrder;
    NumericVector newCase;
    // extend dimension for each case
    for (i = 0; i < cases.size(); i++)
    {
        // compute remaining order
        if (orders[i] + coordinateMaxOrder <= maxOrder)
        {
            remainingOrder = coordinateMaxOrder;
        }
        else
        {
            remainingOrder = maxOrder - orders[i];
        }
        // attach additional dimension
        for (j = 0; j <= remainingOrder; j++)
        {
            newCase = clone(cases[i]);
            newCase.push_back(j);
            newOrders.push_back(orders[i] + j);
            newCases.push_back(newCase);
        }
    }
    // update members
    orders = newOrders;
    cases = newCases;

    // return instance, which automatically converts to reference
    return *this;
}

NumericMatrix VectorEnumerator::getAllCases()
{
    // initialize all cases table
    NumericMatrix allCases(cases.size(), cases[0].size());
    // fill in the cases table
    for (unsigned int i = 0; i < cases.size(); i++)
    {
        allCases(i, _) = cases[i];
    }
    // return all cases
    return allCases;
}

VectorCounter::VectorCounter(NumericVector max, int order)
{
    // define iterators and auxiliary variables
    int i, maxFirstOrder;
    // enumerate all cases
    if (max(0) <= order)
    {
        maxFirstOrder = max(0);
    }
    else
    {
        maxFirstOrder = order;
    }
    VectorEnumerator enumerator(maxFirstOrder);
    for (i = 1; i < max.size(); i++)
    {
        enumerator.extendDimension(max(i), order);
    }
    // get all cases table
    table = enumerator.getAllCases();
    // initialize count and max count
    count = 0;
    maxcount = table.nrow();
}

VectorCounter::VectorCounter(NumericVector max) : VectorCounter(max, sum(max))
{
    //
}

int VectorCounter::getCount()
{
    return count;
}

int VectorCounter::getMaxcount()
{
    return maxcount;
}

NumericMatrix VectorCounter::getTable()
{
    return table;
}

NumericVector VectorCounter::degree()
{
    return table(count, _);
}

void VectorCounter::increase()
{
    count = count + 1;
}

bool VectorCounter::inRange()
{
    if (count < maxcount)
    { // strict inequality since count is zero-base
        return true;
    }
    else
    {
        return false;
    }
}
