#ifndef OPTPOLY_H
#define OPTPOLY_H

#include <Rcpp.h>
using namespace Rcpp;

int choose_cpp(int n, int k);

NumericVector getPrimes_cpp(int n);

NumericMatrix matrixToSparseTriplet_cpp(NumericMatrix mat, bool lowerTriangular);

NumericMatrix vecToMatrix_cpp(NumericVector vectorOfEntries, int dimMatrix, bool lowerTriangular);

NumericMatrix columnEchelon_cpp(NumericMatrix mat);

int degreeToMonomial_cpp(NumericVector degree, NumericVector monomialPrimes);

NumericMatrix monomialsToDegrees_cpp(NumericVector monomials, NumericVector monomialPrimes);

int determineMonomialPosition_cpp(NumericVector degree, NumericVector monomialPrimes, NumericVector momentVector);

double evaluatePolynomial_cpp(NumericVector point, NumericVector coefs, NumericMatrix degrees);

List computeDerivative_cpp(NumericVector coefs, NumericMatrix degrees, int dim, NumericVector monomialPrimes);

List createMosekSdpCoefficientMatrixFromDegrees_cpp(NumericVector coefs, NumericMatrix degrees, NumericVector monomialPrimes, NumericMatrix momentMatrixSparse);

List createMosekSdpCoefficientMatrixFromMonomials_cpp(NumericVector coefs, NumericVector monomials, NumericMatrix momentMatrixSparse);

List createMosekSdpModelSkeletonLasserre_cpp(NumericVector nvar, NumericVector order, NumericVector isConstrained, NumericVector isRectangular, 
                                             NumericMatrix bounds, NumericVector radius, NumericVector monomialPrimes, NumericMatrix momentMatrixSparse);

List createMosekSdpModelSkeletonNieetal_cpp(NumericVector nvar, NumericVector orderObj, NumericVector orderMom,
                                            List gradientObj, NumericVector monomialPrimes, NumericMatrix momentMatrixSparse);

List createMosekQuadraticModelSkeleton_cpp(NumericVector coefs, NumericMatrix degrees, NumericVector monomialPrimes);

#endif