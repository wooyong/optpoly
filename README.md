# optpoly - global polynomial optimization

### Introduction

This package solves global polynomial optimization problems using semidefinite programming (SDP) approaches proposed by Lasserre (2001) and Nie et al (2006).

### Installation

This package uses [Mosek](https://www.mosek.com/)'s R interface `Rmosek` to solve the SDPs, which must be installed before using `optpoly`.

Mosek offers free license to academic users. To install Mosek and its R interface, check the following [installation guide](https://docs.mosek.com/9.0/rmosek/install-interface.html).

To install **optpoly**, type

```r
install.packages("https://wooyong.github.io/data/packages/optpoly_1.1.2.tar.gz", repos=NULL, type="source")
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
sol = optpoly("min", coefs, degrees, opt=NULL)
```

`sol` has the following arguments:

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
 [1,]  1.000000e+00 -9.585283e-18  7.616696e-17  8.071709e-03 -6.402697e-02
 [2,] -9.585283e-18  8.071719e-03 -6.402697e-02 -8.232172e-20  6.139599e-19
 [3,]  7.616696e-17 -6.402697e-02  5.078843e-01  6.531797e-19 -4.879001e-18
 [4,]  8.071709e-03 -8.232172e-20  6.531797e-19  6.547742e-05 -5.169117e-04
 [5,] -6.402697e-02  6.139599e-19 -4.879001e-18 -5.169117e-04  4.099519e-03
 [6,]  5.078843e-01 -4.869281e-18  3.869230e-17  4.099509e-03 -3.251829e-02
 [7,] -5.368288e-20  6.546717e-05 -5.169117e-04 -8.760306e-22  3.334345e-21
 [8,]  6.132222e-19 -5.169117e-04  4.099509e-03  5.922014e-21 -3.945749e-20
 [9,] -8.596957e-18  4.099509e-03 -3.251829e-02 -7.561097e-20  5.562050e-19
[10,]  4.168205e-17 -3.251829e-02  2.579464e-01  3.637811e-19 -2.679651e-18
               [,6]          [,7]          [,8]          [,9]         [,10]
 [1,]  5.078843e-01 -5.368288e-20  6.132222e-19 -8.596957e-18  4.168205e-17
 [2,] -4.869281e-18  6.546717e-05 -5.169117e-04  4.099509e-03 -3.251829e-02
 [3,]  3.869230e-17 -5.169117e-04  4.099509e-03 -3.251829e-02  2.579464e-01
 [4,]  4.099509e-03 -8.760306e-22  5.922014e-21 -7.561097e-20  3.637811e-19
 [5,] -3.251829e-02  3.334345e-21 -3.945749e-20  5.562050e-19 -2.679651e-18
 [6,]  2.579464e-01 -2.733149e-20  3.115810e-19 -4.365337e-18  2.117568e-17
 [7,] -2.733149e-20  1.639651e-06 -4.565089e-06  5.082460e-05 -2.731624e-04
 [8,]  3.115810e-19 -4.565089e-06  5.083485e-05 -2.731624e-04  8.035985e-03
 [9,] -4.365337e-18  5.082460e-05 -2.731624e-04  8.035995e-03 -1.926994e-02
[10,]  2.117568e-17 -2.731624e-04  8.035985e-03 -1.926994e-02  2.231173e+00

$varDim
[1] 2

$order
[1] 6

$certificate
[1] TRUE

$rank
[1] 2

$hierarchy
[1] 1
```

`objective_primal` and `objective_dual` record values of the SDP primal and dual solutions. They are equal when no error occured in solving the SDP. When `certificate=TRUE`, this value is the exact global optimum.

`sdpstatus` and `solstatus` record status of the SDP solution. These are produced by [Mosek](https://www.mosek.com/).

`certificate` indicates the certificate of optimality. If `TRUE`, then `objective_primal` (which equals to `objective_dual`) is the exact global minimum. If `FALSE`, then `min(objective_primal, objective_dual)` is a lower bound for the global minimum.

`rank` represents the number of optimizers, if `certificate=TRUE`. The optimizers can be extracted by the `extractSolution` function in **optpoly**.

`hierarchy` indicates the number of SDP models solved. This is similar to the number of iterations in local optimization algorithms. By default, `optpoly` solves only **one** SDP. This may result in `certificate=FALSE`, in which case one should increase the number of iterations to obtain exact solution. To increase it, specify, for example, `opt=list(hierarchy=3)`.

For more details, including polynomial optimization on a bounded domain, type `?optpoly` in R.

To extract optimizers from the SDP solution, type

```r
extractSolution(sol)
```

Each row of the function's output matrix represents an optimizer. For example, below tells that `(-0.08984222, 0.7126597)` and `(0.08984222 -0.7126597)` are global minimizers of `f(x1,x2)`.

```r
            [,1]       [,2]
[1,] -0.08984222  0.7126597
[2,]  0.08984222 -0.7126597
```

For more details, type `?extractSolution` in R.

### References

Lasserre, Jean B (2001). Global optimization with polynomials and the problem of moments, *SIAM Journal on optimization*, **11**, 796-817.

Nie et al (2006). Minimizing polynomials via sum of squares over the gradient ideal, *Mathematical programming*, **106**, 587-606.

Lasserre, Jean B (2015). *An introduction to polynomial and semi-algebraic optimization*, Cambridge University Press.
