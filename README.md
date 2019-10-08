# optpoly - global polynomial optimization

### Introduction

This package solves global optimization problems of polynomials. It uses semidefinite programming (SDP) approach for global polynomial optimization proposed in Lasserre (2001) and Nie et al (2006).

Given a polynomial, it first creates an R object that describes the SDP model. Then the object is supplied to an SDP solver for the result.

### Installation

This package uses [Mosek](https://www.mosek.com/)'s R interface to solve SDP models, which should be installed before use. If one prefers other SDP solvers, one may take the SDP model object and supply it to the SDP solver of one's choice.

To install Mosek and its R interface, check the following [installation guide](https://docs.mosek.com/9.0/rmosek/install-interface.html).

To install **optpoly**, type

```r
install.packages("https://wooyong.github.io/data/packages/optpoly_2019.10.08.tar.gz", repos=NULL, type="source")
```

Alternatively, to install **optpoly** directly from source on Github, type

```r
# if devtools package is not installed, install it
install.packages("devtools")

# install optpoly from Github using devtools command
library(devtools)
install_github("wooyong/optpoly")
```

To load the package, simply type

```r
library(optpoly)
```

### Usage

Polynomials in **optpoly** are represented by degrees of monomials and their coefficients.

For example, `f(x1, x2) = 4x1^2 - 2.1x1^4 + (1/3)x1^6 + x1x2 - 4x2^2 + 4x2^4` is represented by

```r
degrees = matrix(c(2,0,  # x1^2
                   4,0,  # x1^4
                   6,0,  # x1^6
                   1,1,  # x1x2
                   0,2,  # x2^2
                   0,4), # x2^4
                   nrow=6, ncol=2, byrow=TRUE)
coefs   = c(4, -2.1, 1/3, 1, -4, 4)
```

To globally minimize this polynomial, type

```r
require(optpoly)
optpoly("min", coefs, degrees)
```

The output is as follows.

```r
$objective_primal
[1] -1.031628

$objective_dual
[1] -1.031628

$sdpstatus
[1] "MSK_RES_OK: No error occurred."

$solstatus
[1] "OPTIMAL"

$momentmatrix
               [,1]          [,2]          [,3]          [,4]          [,5]
 [1,]  1.000000e+00  6.444656e-17 -5.053775e-16  8.079834e-03 -6.404979e-02
 [2,]  6.444656e-17  8.079879e-03 -6.404979e-02  3.245976e-19 -4.128110e-18
 [3,] -5.053775e-16 -6.404979e-02  5.077902e-01 -2.431774e-18  3.227131e-17
 [4,]  8.079834e-03  3.245976e-19 -2.431774e-18  6.881513e-05 -5.187924e-04
 [5,] -6.404979e-02 -4.128110e-18  3.227131e-17 -5.187924e-04  4.103023e-03
 [6,]  5.077902e-01  3.282249e-17 -2.573479e-16  4.102977e-03 -3.252392e-02
 [7,]  7.191354e-19  6.876954e-05 -5.187924e-04  3.826438e-20 -8.657027e-20
 [8,] -4.259957e-18 -5.187924e-04  4.102977e-03 -4.938537e-20  2.918671e-19
 [9,]  3.284430e-17  4.102977e-03 -3.252392e-02  1.682414e-19 -2.106020e-18
[10,] -2.575008e-16 -3.252392e-02  2.578509e-01 -1.243436e-18  1.644403e-17
               [,6]          [,7]          [,8]          [,9]         [,10]
 [1,]  5.077902e-01  7.191354e-19 -4.259957e-18  3.284430e-17 -2.575008e-16
 [2,]  3.282249e-17  6.876954e-05 -5.187924e-04  4.102977e-03 -3.252392e-02
 [3,] -2.573479e-16 -5.187924e-04  4.102977e-03 -3.252392e-02  2.578509e-01
 [4,]  4.102977e-03  3.826438e-20 -4.938537e-20  1.682414e-19 -1.243436e-18
 [5,] -3.252392e-02 -8.657027e-20  2.918671e-19 -2.106020e-18  1.644403e-17
 [6,]  2.578510e-01  3.835813e-19 -2.177101e-18  1.672798e-17 -1.311232e-16
 [7,]  3.835813e-19  1.182823e-05 -8.378655e-06  3.523024e-05 -2.636468e-04
 [8,] -2.177101e-18 -8.378655e-06  3.527582e-05 -2.636468e-04  2.083519e-03
 [9,]  1.672798e-17  3.523024e-05 -2.636468e-04  2.083565e-03 -1.651536e-02
[10,] -1.311232e-16 -2.636468e-04  2.083519e-03 -1.651536e-02  1.309342e-01

$varDim
[1] 2

$order
[1] 6

$largeRank
[1] 2

$smallRank
[1] 1

$certificate
[1] FALSE

$hierarchy
[1] 1
```

`objective_primal` and `objective_dual` record global minimum (they are equal when there is no error). `sdpstatus` and `solstatus` record status of SDP solution, which is produced by [Mosek](https://www.mosek.com/). 


### References

Lasserre, Jean B (2001). Global optimization with polynomials and the problem of moments, *SIAM Journal on optimization*, **11**, 796-817.

Nie et al (2006). Minimizing polynomials via sum of squares over the gradient ideal, *Mathematical programming*, **106**, 587-606.

Lasserre, Jean B (2015). *An introduction to polynomial and semi-algebraic optimization*, Cambridge University Press.