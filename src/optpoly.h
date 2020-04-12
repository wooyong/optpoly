#ifndef OPTPOLY_H
#define OPTPOLY_H

#include <Rcpp.h>
using namespace Rcpp;

int choose_cpp(int n, int k);

NumericVector getPrimes_cpp(int n);

NumericMatrix matrixToSparseTriplet_cpp(const NumericMatrix& mat, bool lowerTriangular);

NumericMatrix vecToMatrix_cpp(const NumericVector& vectorOfEntries, int dimMatrix, bool lowerTriangular);

NumericMatrix columnEchelon_cpp(const NumericMatrix& mat);

int degreeToMonomial_cpp(const NumericVector& degree, const NumericVector& monomialPrimes);

NumericMatrix monomialsToDegrees_cpp(const NumericVector& monomials, const NumericVector& monomialPrimes);

int determineMonomialPosition_cpp(const NumericVector& degree, const NumericVector& monomialPrimes, const NumericVector& momentVector);

double evaluatePolynomial_cpp(const NumericVector& point, const NumericVector& coefs, const NumericMatrix& degrees);

List computeDerivative_cpp(const NumericVector& coefs, const NumericMatrix& degrees, int dim, const NumericVector& monomialPrimes);

List createMosekSdpCoefficientMatrixFromDegrees_cpp(NumericVector coefs, NumericMatrix degrees, NumericVector monomialPrimes, NumericMatrix momentMatrixSparse);

List createMosekSdpCoefficientMatrixFromMonomials_cpp(NumericVector coefs, NumericVector monomials, NumericMatrix momentMatrixSparse);

List createMosekSdpModelSkeletonLasserre_cpp(NumericVector nvar, NumericVector order, NumericVector isConstrained, NumericVector isRectangular, 
                                             NumericMatrix bounds, NumericVector radius, NumericVector monomialPrimes, NumericMatrix momentMatrixSparse);

List createMosekSdpModelSkeletonNieetal_cpp(NumericVector nvar, NumericVector orderObj, NumericVector orderMom,
                                            List gradientObj, NumericVector monomialPrimes, NumericMatrix momentMatrixSparse);

List createMosekQuadraticModelSkeleton_cpp(NumericVector coefs, NumericMatrix degrees, NumericVector monomialPrimes);

List createGurobiQuadraticModelSkeleton_cpp(NumericVector coefs, NumericMatrix degrees, NumericVector monomialPrimes);

#endif