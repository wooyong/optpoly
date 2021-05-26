# optpoly - global polynomial optimization

### Introduction

This package solves global polynomial optimization problems using semidefinite programming (SDP) approach proposed by Lasserre (2001) and Nie et al (2006).

### Installation

This package uses [Mosek](https://www.mosek.com/) (version 8 or higher) and its R interface `Rmosek` to solve the SDPs, which must be installed before using `optpoly`.

Mosek offers free academic license. To install Mosek and its R interface, check the following [installation guide](https://docs.mosek.com/9.2/rmosek/install-interface.html).

To install **optpoly**, download the package from [Github release tab](https://github.com/wooyong/optpoly/releases), or copy-and-paste the following to your R console:

```r
install.packages("https://github.com/wooyong/optpoly/releases/download/untagged-6e7c941cd863955bc632/optpoly_1.3.0.tar.gz", repos=NULL, type="source")
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

### Polynomial Optimization

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

`objective_primal` and `objective_dual` are primal and dual solutions of the SDP. They are equal if the SDP algorithm has converged. This equal value is the exact global optimum if `certificate=TRUE`.

`sdpstatus` and `solstatus` are produced by [Mosek](https://www.mosek.com/). They record status of the SDP solution.

`certificate` is the **certificate of optimality**. If `TRUE`, then `objective_primal` (which equals `objective_dual`) is the exact global minimum. If `FALSE`, then `min(objective_primal, objective_dual)` is a lower bound for the global minimum.

`rank` is the number of optimizers, provided `certificate=TRUE`. Optimizers can be extracted by the `extractSolution` function explained in the next section.

`hierarchy` is the number of SDPs solved to reach the current output. This is analogous to the number of iterations in other optimization algorithms. By default, `optpoly` solves only **one** SDP. This may result in `certificate=FALSE`, in which case one should increase it until one gets `certificate=TRUE`. To increase it, specify e.g. `opt=list(hierarchy=3)` in the option.

For more details, including polynomial optimization on a bounded domain, type `?optpoly` in R.

### Extracting Solutions

To extract optimizers from the SDP solution, type

```r
extractSolution(sol)

            [,1]       [,2]
[1,] -0.08984222  0.7126597
[2,]  0.08984222 -0.7126597

```

Each row of the output matrix represents an optimizer. The above tells that `(-0.08984222, 0.7126597)` and `(0.08984222 -0.7126597)` are global minimizers of `f(x1,x2)`.

For more details, type `?extractSolution` in R.

### System of Polynomial Equations

The SDP approach of polynomial optimization can also be used to solve a system of polynomial equations. The SDP approach directly searchs for real-valued solutions, which is computationally cheaper than algorithms that search for both real and complex solutions.

To solve a system of polynomial equations, the SDP approach minimizes an auxiliary function (e.g. a zero function) with the system of equations as constraints. It then extract solutions from it.

This is implemented by `solvepoly` function in **optpoly**. As an example, to solve a system `x^2+y^2-20=0` and `x+y-6=0`, type

```r
# define system of equations
eqList = list(
    # x^2 + y^2 - 20 = 0
    list(degrees = matrix(c(2,0,  # x^2
                            0,2,  # y^2
                            0,0), # constant
                            nrow=3, ncol=2, byrow=TRUE),
         coefs   = c(1, 1, -20)),
    # x + y - 6 = 0
    list(degrees = matrix(c(1,0,
                            0,1,
                            0,0),
                            nrow=3, ncol=2, byrow=TRUE),
        coefs    = c(1, 1, -6))
)

# solve system of equations
require(optpoly)
sol = solvepoly(eqList, opt=list(hierarchy=2))
extractSolution(sol)

     [,1] [,2]
[1,]    2    4
[2,]    4    2
```

`sol` has the same output structure as `optpoly`. It is an intermediate output from minimizing a zero function with the system of equations as constraints. `sol` is then taken to `extractSolution` function to retrieve the solutions of the system of polynomial equations.

### References

Lasserre, Jean B (2001). Global optimization with polynomials and the problem of moments, *SIAM Journal on optimization*, **11**, 796-817.

Nie et al (2006). Minimizing polynomials via sum of squares over the gradient ideal, *Mathematical programming*, **106**, 587-606.

Lasserre, Jean B (2015). *An introduction to polynomial and semi-algebraic optimization*, Cambridge University Press.
