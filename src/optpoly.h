#ifndef OPTPOLY_H
#define OPTPOLY_H

#include <Rcpp.h>
using namespace Rcpp;

int choose_cpp(int n, int k);

int degreeToMonomial_cpp(NumericVector degree, NumericVector monomialPrimes);

NumericMatrix monomialsToDegrees_cpp(NumericVector monomials, NumericVector monomialPrimes);

int determineMonomialPosition_cpp(NumericVector degree, NumericVector monomialPrimes, NumericVector momentVector);

double evaluatePolynomial_cpp(NumericVector point, NumericVector coefs, NumericMatrix degrees);

List computeDerivative_cpp(NumericVector coefs, NumericMatrix degrees, int dim, NumericVector monomialPrimes);

List createMosekSdpCoefficientMatrixFromDegrees_cpp(NumericVector coefs, NumericMatrix degrees, NumericVector monomialPrimes, NumericMatrix momentMatrixSparse);

List createMosekSdpCoefficientMatrixFromMonomials_cpp(NumericVector coefs, NumericVector monomials, NumericMatrix momentMatrixSparse);

List createMosekSdpModelSkeleton_cpp(NumericVector nvar, NumericVector order, NumericVector isConstrained, NumericVector radius,
                                     NumericVector monomialPrimes, NumericMatrix momentMatrixSparse);

List createMosekSdpModelSkeletonWithGradientIdeals_cpp(NumericVector nvar, NumericVector orderObj, NumericVector orderMom,
                                                       List gradientObj, NumericVector monomialPrimes, NumericMatrix momentMatrixSparse);

#endif