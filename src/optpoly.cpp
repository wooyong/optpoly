#include "optpoly.h"
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
int choose_cpp(int n, int k) {

    /*

    Computes n choose k

    */

    int result = 1;
    int j;

    for(j=0; j<k; j++) {
        result = result * (n-j) / (j+1);
    }

    return result;
}

// [[Rcpp::export]]
int degreeToMonomial_cpp(NumericVector degree, NumericVector monomialPrimes) {

    // initialize identifier
    double monomialID;
    monomialID = 1;

    // compute identifier
    int k;
    for(k=1; k<=monomialPrimes.size(); k++) {
        monomialID = monomialID * pow(monomialPrimes(k-1), degree(k-1));
    }

    return static_cast<int>(monomialID);
}

// [[Rcpp::export]]
NumericMatrix monomialsToDegrees_cpp(NumericVector monomials, NumericVector monomialPrimes) {

    // read dimension
    int P = monomials.size();
    int L = monomialPrimes.size();

    // initialize degrees matrix
    NumericMatrix degrees(P, L);

    // initialize auxiliary variable
    int monomial, factor;

    // perform factorization to extract degrees
    for(int j=0; j<P; j++) {
        monomial = monomials(j);
        for(int k=0; k<L; k++) {
            factor = monomialPrimes(k);
            while(monomial % factor == 0) {
                monomial = monomial / factor;
                degrees(j,k) = degrees(j,k) + 1;
            }
        }
    }

    // return degrees matrix
    return degrees;
}

// [[Rcpp::export]]
int determineMonomialPosition_cpp(NumericVector degree, NumericVector monomialPrimes, NumericVector momentVector) {

    /*

    Finds the location of the monomial represented by _degree_ in _momentVector_.

    The index begins with 0. That is, if _degree_ = [0, ..., 0] and the first entry of _momentVector_ represents the constant,
        then this function returns 0.

    It first computes the number representing the monomial from the list of prime numbers in _monomialPrimes_
        and then find its position in _momentVector_.

    *** not used for now, kept it for future use ***

    */

    // get monomial ID
    int monomialID = degreeToMonomial_cpp(degree, monomialPrimes);

    // return vector location. Subtract 1 to reflect 0-base.
    NumericVector position(1);
    position = match(NumericVector::create(monomialID), momentVector) - 1;
    return position(0);
}

// [[Rcpp::export]]
double evaluatePolynomial_cpp(NumericVector point, NumericVector coefs, NumericMatrix degrees) {

    // read dimension of variables
    int dim = degrees.ncol();

    // read number of monomials
    int P = degrees.nrow();

    // copy coefficient vector
    NumericVector value = clone(coefs);

    // multiply monomial value to coefficients
    for(int k=0; k<dim; k++) {
        for(int j=0; j<P; j++) {
            value(j) = value(j) * pow(point(k), degrees(j,k));
        }
    }

    // return sum of monomials
    return sum(value);
}

// [[Rcpp::export]]
List computeDerivative_cpp(NumericVector coefs, NumericMatrix degrees, int dim, NumericVector monomialPrimes) {

    // read number of monomials
    int P = degrees.nrow();

    // initialize derivative
    NumericVector dCoefs   = clone(coefs);
    NumericMatrix dDegrees = clone(degrees);

    // differentiate each monomial
    for(int j=0; j<P; j++) {
        if(dDegrees(j,dim-1) > 0) {
            dCoefs(j) = dCoefs(j) * dDegrees(j,dim-1);
            dDegrees(j,dim-1) = dDegrees(j,dim-1) - 1;
        } else {
            dCoefs(j) = 0;
        }
    }

    // return derivative
    return List::create(_["coefs"] = dCoefs,
                        _["degrees"] = dDegrees);
}

// [[Rcpp::export]]
List createMosekSdpCoefficientMatrixFromDegrees_cpp(NumericVector coefs, NumericMatrix degrees, NumericVector monomialPrimes, NumericMatrix momentMatrixSparse) {

    /*

    Creates the polynomial coefficient matrix in sparse triplet form.

    The return value is the input of the SDP model as the coefficients of the objective function.

    Since we intends to use it for Mosek and Mosek copies the lower triangular part to the upper triangular part automatically,
        we reduce the coefficient of the off-diagonals by half.

    */

    // read dimensions
    int P = degrees.nrow();
    int L = monomialPrimes.size();

    // create monomial vector
    NumericVector monomialVector(P);

    // define iterators
    int k, m;
    double monomialID;

    // fill in the monomial vector
    for(k=0; k<P; k++) {
        monomialID = 1;
        for(m=0; m<L; m++) {
            monomialID = monomialID * pow(monomialPrimes(m), degrees(k,m));
        }
        monomialVector(k) = monomialID;
    }

    // compute and return coefficient matrix from monomial vector
    return createMosekSdpCoefficientMatrixFromMonomials_cpp(coefs, monomialVector, momentMatrixSparse);
}

// [[Rcpp::export]]
List createMosekSdpCoefficientMatrixFromMonomials_cpp(NumericVector coefs, NumericVector monomials, NumericMatrix momentMatrixSparse) {

    /*

    Creates the polynomial coefficient matrix in sparse triplet form.

    The return value is the input of the SDP model as the coefficients of the objective function.

    Since we intends to use it for Mosek and Mosek copies the lower triangular part to the upper triangular part automatically,
        we reduce the coefficient of the off-diagonals by half.

    */

    // read dimensions
    int P = monomials.size();
    int J = momentMatrixSparse.nrow();

    // initialize coefficient vector
    NumericVector coefVector(J);

    // define iterators and auxiliary variables
    int k;
    NumericVector loc(1);

    // fill in the coefficient matrix
    for(k=0; k<P; k++) {
        // subtract 1 to reflect 0-base
        loc = match(NumericVector::create(monomials(k)), momentMatrixSparse(_,2)) - 1;
        // if diagonal, add coefficient. If off-diagonal, add half of the coefficient
        if(momentMatrixSparse(loc(0),0) == momentMatrixSparse(loc(0),1)) {
            coefVector(loc(0)) += coefs(k);
        } else {
            coefVector(loc(0)) += 0.5 * coefs(k);
        }
    }

    // return coefficient matrix
    return List::create(_["j"] = rep(1, J),
                        _["k"] = momentMatrixSparse(_,0),
                        _["l"] = momentMatrixSparse(_,1),
                        _["v"] = coefVector);
}

// [[Rcpp::export]]
List createMosekSdpModelSkeleton_cpp(NumericVector nvar, NumericVector order, NumericVector isConstrained, NumericVector radius,
                                     NumericVector monomialPrimes, NumericMatrix momentMatrixSparse) {

    /*

    The SDP problem in the standard form is defined by

        minimize        sum_{j=1}^{n} c_j * x_j + sum_{j=1}^{p} tr<bar C_j, bar X_j> + c^f

        by choosing     x_j,        j=1,...,n,      the real numbers,

                        bar X_j,    j=1,...,p,      the r_j * r_j symmetric matrices,

        subject to      l_i^c <= sum_{j=1}^{n} a_{ij} * x_j + sum_{j=1}^{p} tr<bar A_{ij}, bar X_j> <= u_i^c,   i = 1,...,m,

                        l_j^x <= x_j <= u_j^x,                                                                  j = 1,...,n,

                        bar X_j is positive semidefinite, j = 1,...,p.

    In Mosek's R interface, the model is specified using a List containing:

        CharacterVector sense   = "min" or "max",
        NumericVector   c       = (c_1,...,c_n),
        NumericMatrix   A       = [a_{ij}],   i=1,...,m,   j=1,...,n,
        NumericMatrix   bc      = rbind(blc = (l_1^c,...,l_m^c),
                                        buc = (u_1^c,...,u_m^c)),
        NumericMatrix   bx      = rbind(blx = (l_1^x,...,l_n^x),
                                        bux = (u_1^x,...,u_n^x)),
        NumericVector   bardim  = (r_1,...,r_n),
        List            barc    = list(j=(1,...,1,                         2,...,2,                         ..., p,...,p),
                                       k=(rows.in.SparseMatrix(bar C_1),   rows.in.SparseMatrix(bar C_2),   ..., rows.in.SparseMatrix(bar C_p)),
                                       l=(cols.in.SparseMatrix(bar C_1),   cols.in.SparseMatrix(bar C_2),   ..., cols.in.SparseMatrix(bar C_p)),
                                       v=(values.in.SparseMatrix(bar C_1), values.in.SparseMatrix(bar C_2), ..., values.in.SparseMatrix(bar C_p))),
        List            barA    = list(i=(1,...,1,                                                                      ..., m,...,m),
                                       j=(1,...,1,                             ..., p,...,p,                            ..., 1,...,1,                            ..., p,...,p)
                                       k=(rows.in.SparseMatrix(bar A_{11}),    ..., rows.in.SparseMatrix(bar A_{1p}),   ..., rows.in.SparseMatrix(bar A_{m1}),   ..., rows.in.SparseMatrix(bar A_{mp})),
                                       l=(cols.in.SparseMatrix(bar A_{11}),    ..., cols.in.SparseMatrix(bar A_{1p}),   ..., cols.in.SparseMatrix(bar A_{m1}),   ..., cols.in.SparseMatrix(bar A_{mp})),
                                       v=(values.in.SparseMatrix(bar A_{11})), ..., values.in.SparseMatrix(bar A_{1p}), ..., values.in.SparseMatrix(bar A_{m1}), ..., values.in.SparseMatrix(bar A_{mp})))

        where   rows.in.SparseMatrix is the vector of the row indices in the sparse matrix's triplet representation,
                cols.in.SparseMatrix is the vector of the column indices in the sparse matrix's triplet representation,
                values.in.SparseMatrix is the vector of the value indices in the sparse matrix's triplet representation.

    */

    /*

    The polynomial optimization problem in SDP is defined by the following.

    Suppose we want to optimize a 4th order polynomial in (x,y).

    Let monomialVector = (1 x y x2 xy y2)' be the column vector of all monomials in (x,y) up to 2nd order,
        enumerated from the lowest order (i.e. zero order) to the highest order.

    Let M = outer.product(monomialVector, monomialVector), which yields

        M = [1    x    y    x2   xy   y2
             x    x2   xy   x3   x2y  xy2
             y    xy   y2   x2y  xy2  y3
             x2   x3   x2y  x4   x3y  x2y2
             xy   x2y  xy2  x3y  x2y2 xy3
             y2   xy2  y3   x2y2 xy3  y4  ].

    Let C be the matrix of coefficients corresponding to each entry of M.
        In principle, the C matrix doesn't need to be symmetric, 
        e.g. one may record coefficients only using the lower-triangular proportion of M.

    However, Mosek assumes C is symmetric and copies the lower triangular part of C to its the upper triangular part.
        Therefore, one must reduce its coefficient by half when off-diagonal.

    The SDP formulation of the polynomial optimization is done by
        selecting the value of the monomials arranged in the matrix form according to M.

    The SDP problem for polynomial optimization is defined by

        minimize        tr<C, M>

        by choosing     M,

        subject to      M is positive semidefinite.

    In the standard form (X in the below is M in the above):

        n = 1, p = 1, r_1 = 6,

        minimize        0 * x_1 + tr<C, X>

        by choosing     x_1 and X,

        subject to      X[1,1] = 1,             The value of the zero-order monomial must be 1.
                        X[4,1] - X[2,2] = 0,    The monomial 'x2' appears multiple times in M, and the monomial values must be equal.
                        X[5,1] - X[3,2] = 0,    The monomial 'xy' appears multiple times in M, and the monomial values must be equal.
                        ...,
                        X[6,2] - X[5,3] = 0,    The monomial 'xy2' appears multiple times in M, and the monomial values must be equal.
                        X is positive semidefinite.

    Mosek identifies the matrix variable X by its lower triangular part,
        so the equivalence constraints of the monomial values are listed only for the lower triangular part.

    Higher order SDPs are introduced by using the monomialVector of higher order.
        In the unconstrained case, this does not make any difference.
        In fact, it will be worse because the numerical complexity increases and the SDP algorithms may not converge.
        However, in the constrained case, it will have an effect: the SDP solution will be a better approximation to the true solution.

    */

    /*

    Suppose we want to perform the above polynomial optimization on a bounded domain. Then we need to introduce additional constraints.
        Suppose the domain is {(x,y) | x^2 + y^2 <= R^2} for some positive real number R.

    Let M' be the upper left part of M where the entries are up to 2nd order, i.e.

        M' = [1    x    y
              x    x2   xy
              y    xy   y2].

    Generally, if M has up to A-th order and the constraint is of B-th order, then we choose M' so that it contains up to (A-B)-th order.
        The M' matrix is called the localizing matrix.

    Define Y = R^2 * M' - x^2 * M' - y^2 * M', that is,

        Y[1,1] = R^2 - x^2 - y^2,
        Y[2,1] = R^2*x - x^3 - x*y^2,
        ...,
        Y[3,3] = R^2*y^2 - x^2*y^2 - y^4.

    Then the additional constraint is that the matrix Y, which is the matrix of linear combinations of monomial values, is positive semidefinite.

    So, in the standard form, the SDP problem for polynomial optimization on the bounded domain is

        n = 1, p = 2, r_1 = 6, r_2 = 3,

        minimize        0 * x_1 + tr<C, X1> + tr<0, X2>

        by choosing     x_1, X1, and X2

        subject to      X1[1,1] = 1,                                            The value of the zero-order monomial must be 1.
                        X1[4,1] - X1[2,2] = 0,                                  The monomial 'x2' appears multiple times in M, and their values must be equal.
                        X1[5,1] - X1[3,2] = 0,                                  The monomial 'xy' appears multiple times in M, and their values must be equal.
                        ...,
                        X1[6,2] - X1[5,3] = 0,                                  The monomial 'xy2' appears multiple times in M, and their values must be equal.
                        X2[1,1] - R^2 * X1[1,1] + X1[4,1] + X1[6,1] = 0,        Y[1,1] := R^2 * 1 - x^2 - y^2.
                        X2[2,1] - R^2 * X1[2,1] + X1[4,2] + X1[5,3] = 0,        Y[2,1] := R^2 * x - x^3 - x*y^2.
                        ...,
                        X2[3,3] - R^2 * X1[3,3] + X1[6,4] + X1[6,6] = 0,        Y[3,3] := R^2 * y^2 - x^2*y^2 - y^4.
                        X1 is positive semidefinite,
                        X2 is positive semidefinite.

    */

    /*

    Now we write the SDP formulation of the polynomial optimization problem
        in the structure that Mosek can understand.

    */

    /*

    We first count the number of constraints.

    For brevity of notation in the exposition, let
        dim  = dimMomentMatrix,
        dim' = dimLocalizingMatrix.

    First, there is a constraint that X[1,1] = 1.

    In addition, there are
        dim + choose(dim,2) = dim*(dim+1)/2     lower triangular entries in momentMatrix, and
        choose(nvar+order, nvar)                unique monomials.

    So there are 
        dim*(dim+1)/2 - choose(nvar+order, nvar)
    number of constraints that impose equivalence of monomial values.

    Also, if the problem is on a bounded domain,
        we additionally have dim' + choose(dim',2) = dim'*(dim'+1)/2 constraints.

    To code the constraints in the structure that Mosek can understand, 
        we also need to compute the number of nonzero entries of the coefficient matrix
        appearing in the constraints.

    The number of nonzero entries are computed as follows:
        First, there is one nonzero entry in the constraint which is X[1,1] = 1.
        Second, there are two nonzero entries for each equivalence-of-monomial-values constraints.
        Third, there are (2+nvar) nonzero entries for each localizing matrix constraints.

    */

    // record nearest even degree
    int evenOrder = 2 * floor( (order(0) + 1) / 2 );

    // record if the problem is on a bounded domain
    int iConstrained = isConstrained(0);

    // compute the number of monomials up to _order_ -th order in _nvar_ variables
    int nMonomials = choose_cpp(static_cast<int>(nvar(0) + evenOrder), static_cast<int>(nvar(0)));

    // compute the dimension of momentMatrix. Note that it equals to the number of monomials up to as.integer( (_order_+1)/2 ).
    int dimMomentMatrix = choose_cpp(static_cast<int>(nvar(0) + evenOrder/2), static_cast<int>(nvar(0)));
    int dimLocalizingMatrix = choose_cpp(static_cast<int>(nvar(0) + evenOrder/2 - 1), static_cast<int>(nvar(0)));

    // number of constraints
    int entriesInMomentMatrix = dimMomentMatrix * (dimMomentMatrix + 1) / 2;
    int entriesInLocalizingMatrix = dimLocalizingMatrix * (dimLocalizingMatrix + 1) / 2;
    int nConstraints = 1 + 
                       entriesInMomentMatrix - nMonomials + 
                       iConstrained * entriesInLocalizingMatrix;

    // number of nonzero coefficient matrix entries in the constraints
    int nNonzeroEntries = 1 + 
                          2 * ( entriesInMomentMatrix - nMonomials ) +
                          iConstrained * (2 + nvar(0)) * entriesInLocalizingMatrix;

    // dimension of the scalar choice variable (x_1 in the standard form) is 1, so we don't define it.
    ///// int nScalarVariable = 1; /////

    /*

    Now we construct Mosek's SDP program structure.

    We first define each piece of the structure and then combine them as a list at the end.

    */

    // optimization sense, either "max" or "min". Left blank in this function
    CharacterVector sense("assign max or min, depending on the problem.");
    // coefficient on the scalar variable in the objective, which is set to zero by default
    NumericVector c(1);
    // coefficient on the scalar variable in the constraints, which are set to zeros by default
    NumericMatrix A(nConstraints,1);
    // lower and upper bounds for the constraints. The entries in the first column are one (which represents X[1,1] = 1), otherwise zero.
    NumericMatrix bc(2,nConstraints);
    bc(0,0) = 1; bc(1,0) = 1;
    // lower and upper bounds for the scalar variable, which are zeros
    NumericMatrix bx(2,1);

    // dimensions of matrix variables. Number of matrix variables change according to iConstrained = isConstrained(0).
    NumericVector *bardim; // first define pointer and then assign values
    NumericVector bardimUnbounded = NumericVector::create(dimMomentMatrix);
    NumericVector bardimBounded = NumericVector::create(dimMomentMatrix, dimLocalizingMatrix);
    if(iConstrained == 0) {
        bardim = &bardimUnbounded;
    } else {
        bardim = &bardimBounded;
    }

    // coefficients on the matrix variable in the objective.
    // this part is left blank since it is problem-specific and not a "skeleton".
    CharacterVector barc("fill in the coefficients in the sparseMatrix triplet form.");

    /*

    coefficients on the matrix variable in the constraints, barA.

    */

    // we define each entry of barA separately.
    NumericVector barA_i(nNonzeroEntries);
    NumericVector barA_j(nNonzeroEntries);
    NumericVector barA_k(nNonzeroEntries);
    NumericVector barA_l(nNonzeroEntries);
    NumericVector barA_v(nNonzeroEntries);

    // first write the X[1,1] = 1 constraint
    barA_i(0) = 1;
    barA_j(0) = 1;
    barA_k(0) = 1;
    barA_l(0) = 1;
    barA_v(0) = 1;

    // define iterators and counters needed to enumerate the rest of the constraints
    int barAInd = 1;   // position in the bara vectors
    int iInd = 2;      // counter for the constraint number
    int j, k, l;
    NumericVector momentId(1);
    NumericVector momentIdPosition(1);

    /*

    write equivalence of the monomial values

    */

    // enumerate unique elements in momentMatrix
    NumericVector momentVector = momentMatrixSparse(_,2);
    NumericVector uniqueMomentVector = unique(momentVector);

    // impose equivalence for each entry in uniqueMomentVector
    for(k=0; k<uniqueMomentVector.size(); k++) {

        // read moment ID
        momentId = uniqueMomentVector(k);
        // record the earliest position of momentID in momentVector
        momentIdPosition = match(momentId, momentVector) - 1;

        // if moment appears again in the remaining sequence of the momentVector, record the equivalence constraint
        // search for additional appearance if the first appearance is not the last entry of momentVector
        if(momentIdPosition(0) < momentVector.size()) {
            for(l=momentIdPosition(0)+1; l<momentVector.size(); l++) {
                // if there is an additional appearance, record the constraint
                if(momentId(0) == momentVector(l)) {
                    // code the constraint
                    barA_i(barAInd) = iInd;
                    barA_j(barAInd) = 1;
                    barA_k(barAInd) = momentMatrixSparse(momentIdPosition(0),0);
                    barA_l(barAInd) = momentMatrixSparse(momentIdPosition(0),1);
                    // if off-diagonal, deflate the coefficient by half since the same is also copied to the upper-triangular part in Mosek
                    if(barA_k(barAInd) == barA_l(barAInd)) {
                        barA_v(barAInd) = 1.0;
                    } else {
                        barA_v(barAInd) = 0.5;
                    }
                    
                    barA_i(barAInd+1) = iInd;
                    barA_j(barAInd+1) = 1;
                    barA_k(barAInd+1) = momentMatrixSparse(l,0);
                    barA_l(barAInd+1) = momentMatrixSparse(l,1);
                    // if off-diagonal, deflate the coefficient by half since the same is also copied to the upper-triangular part in Mosek
                    if(barA_k(barAInd+1) == barA_l(barAInd+1)) {
                        barA_v(barAInd+1) = -1.0;
                    } else {
                        barA_v(barAInd+1) = -0.5;
                    }
                    iInd++;
                    barAInd = barAInd + 2;
                }
            }
        }
    }

    /*

    If the optimization is on bounded domain, impose localizing matrix constraints.

    We consider the bound that the arguments of the polynomial, denoted by x_1, ..., x_{nvar},
        are within some ball of the origin with radius R:

        x_1^2 + ... + x_{nvar}^2 <= R^2.

    The localizing matrix for this bound is

        R^2 * M' - x_1^2 * M' - ... - x_n^2 * M'.

    */

    // define auxiliary variables used to write the constraints
    NumericVector squaredMomentId(1);
    NumericVector squaredMomentIdPosition(1);
    NumericMatrix squaredMomentIdsPositions(nvar(0), 2);

    if(iConstrained == 1) {
        // write constraint for each entry (k,l) of localizing matrix.
        // to enumerate, we loop over rows of momentMatrixSparse and write the constraint for the localizing matrix part of it.
        for(k=0; k<momentMatrixSparse.nrow(); k++) {
            if(momentMatrixSparse(k,0) <= dimLocalizingMatrix) {
                // search for the position of the term x_j^2 * M' in momentMatrixSparse
                for(j=1; j<=nvar(0); j++) {
                    squaredMomentId = momentMatrixSparse(k,2) * monomialPrimes(j-1) * monomialPrimes(j-1);
                    squaredMomentIdPosition = match(squaredMomentId, momentVector) - 1;
                    squaredMomentIdsPositions(j-1,0) = momentMatrixSparse(squaredMomentIdPosition(0),0);
                    squaredMomentIdsPositions(j-1,1) = momentMatrixSparse(squaredMomentIdPosition(0),1);
                }
                // record constraint
                barA_i(barAInd) = iInd;
                barA_j(barAInd) = 2;
                barA_k(barAInd) = momentMatrixSparse(k,0);
                barA_l(barAInd) = momentMatrixSparse(k,1);
                // if off-diagonal, deflate the coefficient by half since the same is also copied to the upper-triangular part in Mosek
                if(barA_k(barAInd) == barA_l(barAInd)) { barA_v(barAInd) = 1.0; } else { barA_v(barAInd) = 0.5; }
                barA_i(barAInd+1) = iInd;
                barA_j(barAInd+1) = 1;
                barA_k(barAInd+1) = momentMatrixSparse(k,0);
                barA_l(barAInd+1) = momentMatrixSparse(k,1);
                // if off-diagonal, deflate the coefficient by half since the same is also copied to the upper-triangular part in Mosek
                if(barA_k(barAInd+1) == barA_l(barAInd+1)) { barA_v(barAInd+1) = - radius(0) * radius(0); } else { barA_v(barAInd+1) = - radius(0) * radius(0) * 0.5; }
                barAInd = barAInd + 2;
                for(j=1; j<=nvar(0); j++) {
                    barA_i(barAInd) = iInd;
                    barA_j(barAInd) = 1;
                    barA_k(barAInd) = squaredMomentIdsPositions(j-1,0);
                    barA_l(barAInd) = squaredMomentIdsPositions(j-1,1);
                    // if off-diagonal, deflate the coefficient by half since the same is also copied to the upper-triangular part in Mosek
                    if(barA_k(barAInd) == barA_l(barAInd)) { barA_v(barAInd) = 1.0; } else { barA_v(barAInd) = 0.5; }
                    barAInd++;
                }
                iInd++;
            }
        }
    }

    // return the skeleton
    return List::create(_["sense"] = sense,
                        _["c"] = c,
                        _["A"] = A,
                        _["bc"] = bc,
                        _["bx"] = bx,
                        _["bardim"] = *bardim,
                        _["barc"] = barc,
                        _["barA"] = List::create(_["i"] = barA_i,
                                                 _["j"] = barA_j,
                                                 _["k"] = barA_k,
                                                 _["l"] = barA_l,
                                                 _["v"] = barA_v));

}

// [[Rcpp::export]]
List createMosekSdpModelSkeletonWithGradientIdeals_cpp(NumericVector nvar, NumericVector orderObj, NumericVector orderMom,
                                                       List gradientObj, NumericVector monomialPrimes, NumericMatrix momentMatrixSparse) {

    /*

    The SDP problem in the standard form is defined by

        minimize        sum_{j=1}^{n} c_j * x_j + sum_{j=1}^{p} tr<bar C_j, bar X_j> + c^f

        by choosing     x_j,        j=1,...,n,      the real numbers,

                        bar X_j,    j=1,...,p,      the r_j * r_j symmetric matrices,

        subject to      l_i^c <= sum_{j=1}^{n} a_{ij} * x_j + sum_{j=1}^{p} tr<bar A_{ij}, bar X_j> <= u_i^c,   i = 1,...,m,

                        l_j^x <= x_j <= u_j^x,                                                                  j = 1,...,n,

                        bar X_j is positive semidefinite, j = 1,...,p.

    In Mosek's R interface, the model is specified using a List containing:

        CharacterVector sense   = "min" or "max",
        NumericVector   c       = (c_1,...,c_n),
        NumericMatrix   A       = [a_{ij}],   i=1,...,m,   j=1,...,n,
        NumericMatrix   bc      = rbind(blc = (l_1^c,...,l_m^c),
                                        buc = (u_1^c,...,u_m^c)),
        NumericMatrix   bx      = rbind(blx = (l_1^x,...,l_n^x),
                                        bux = (u_1^x,...,u_n^x)),
        NumericVector   bardim  = (r_1,...,r_n),
        List            barc    = list(j=(1,...,1,                         2,...,2,                         ..., p,...,p),
                                       k=(rows.in.SparseMatrix(bar C_1),   rows.in.SparseMatrix(bar C_2),   ..., rows.in.SparseMatrix(bar C_p)),
                                       l=(cols.in.SparseMatrix(bar C_1),   cols.in.SparseMatrix(bar C_2),   ..., cols.in.SparseMatrix(bar C_p)),
                                       v=(values.in.SparseMatrix(bar C_1), values.in.SparseMatrix(bar C_2), ..., values.in.SparseMatrix(bar C_p))),
        List            barA    = list(i=(1,...,1,                                                                      ..., m,...,m),
                                       j=(1,...,1,                             ..., p,...,p,                            ..., 1,...,1,                            ..., p,...,p)
                                       k=(rows.in.SparseMatrix(bar A_{11}),    ..., rows.in.SparseMatrix(bar A_{1p}),   ..., rows.in.SparseMatrix(bar A_{m1}),   ..., rows.in.SparseMatrix(bar A_{mp})),
                                       l=(cols.in.SparseMatrix(bar A_{11}),    ..., cols.in.SparseMatrix(bar A_{1p}),   ..., cols.in.SparseMatrix(bar A_{m1}),   ..., cols.in.SparseMatrix(bar A_{mp})),
                                       v=(values.in.SparseMatrix(bar A_{11})), ..., values.in.SparseMatrix(bar A_{1p}), ..., values.in.SparseMatrix(bar A_{m1}), ..., values.in.SparseMatrix(bar A_{mp})))

        where   rows.in.SparseMatrix is the vector of the row indices in the sparse matrix's triplet representation,
                cols.in.SparseMatrix is the vector of the column indices in the sparse matrix's triplet representation,
                values.in.SparseMatrix is the vector of the value indices in the sparse matrix's triplet representation.

    */

    /*

    The polynomial optimization problem in SDP is defined by the following.

    Suppose we want to optimize a 4th order polynomial in (x,y).

    Let monomialVector = (1 x y x2 xy y2)' be the column vector of all monomials in (x,y) up to 2nd order,
        enumerated from the lowest order (i.e. zero order) to the highest order.

    Let M = outer.product(monomialVector, monomialVector), which yields

        M = [1    x    y    x2   xy   y2
             x    x2   xy   x3   x2y  xy2
             y    xy   y2   x2y  xy2  y3
             x2   x3   x2y  x4   x3y  x2y2
             xy   x2y  xy2  x3y  x2y2 xy3
             y2   xy2  y3   x2y2 xy3  y4  ].

    Let C be the matrix of coefficients corresponding to each entry of M.
        In principle, the C matrix doesn't need to be symmetric, 
        e.g. one may record coefficients only using the lower-triangular proportion of M.

    However, Mosek assumes C is symmetric and copies the lower triangular part of C to its the upper triangular part.
        Therefore, one must reduce its coefficient by half when off-diagonal.

    The SDP formulation of the polynomial optimization is done by
        selecting the value of the monomials arranged in the matrix form according to M.

    The SDP problem for polynomial optimization is defined by

        minimize        tr<C, M>

        by choosing     M,

        subject to      M is positive semidefinite.

    In the standard form (X in the below is M in the above):

        n = 1, p = 1, r_1 = 6,

        minimize        0 * x_1 + tr<C, X>

        by choosing     x_1 and X,

        subject to      X[1,1] = 1,             The value of the zero-order monomial must be 1.
                        X[4,1] - X[2,2] = 0,    The monomial 'x2' appears multiple times in M, and the monomial values must be equal.
                        X[5,1] - X[3,2] = 0,    The monomial 'xy' appears multiple times in M, and the monomial values must be equal.
                        ...,
                        X[6,2] - X[5,3] = 0,    The monomial 'xy2' appears multiple times in M, and the monomial values must be equal.
                        X is positive semidefinite.

    Mosek identifies the matrix variable X by its lower triangular part,
        so the equivalence constraints of the monomial values are listed only for the lower triangular part.

    Higher order SDPs are introduced by using the monomialVector of higher order.
        In the unconstrained case, this does not make any difference.
        In fact, it will be worse because the numerical complexity increases and the SDP algorithms may not converge.
        However, in the constrained case, it will have an effect: the SDP solution will be a better approximation to the true solution.

    */

    /*

    The approach of Nie et al (2006, Math Program Ser A) is to impose the first order conditions via gradient ideals.

    In the above example, let f be the objective polynomial and f_x and f_y be the partial derivative of f with respect to x and y:

        f_x = df / dx, 
        f_y = df / dy.

    Let orderObj be the order of f and orderMom be the max order of the momentMatrix.

    Let M' be the upper left part of M where the entries are up to (orderMom-orderObj)-th order.
        For example, if orderObj = 4 and orderMom = 4 as in the above example, then

        M' = [1].

        If orderObj = 6 and orderMom = 4, then

        M' = [1    x    y
              x    x2   xy
              y    xy   y2].

        With abuse of notation, we also call this a localizing matrix.

    Then the first order conditions imposed via gradient ideals are

        f_x * M' = 0,
        f_y * M' = 0,

        where the right-hand side is the matrix of zeroes.

    For example, f_x * M' = 0 when orderObj=6 means that

        f_x * M'[1,1] = f_x = 0,        the first order condition,
        f_x * M'[1,2] = x * f_x = 0,    the first order condition multiplied by x,
        f_x * M'[1,2] = y * f_x = 0,    the first order condition multiplied by y,
        ...
        f_x * M'[3,3] = y2 * f_x = 0,   the first order condition multiplied by y2.

    We represent this condition in SDP as follows. For example, if f_x = 3x^3 + 2x^2y + y + 1,

        3 * M[4,2] + 2 * M[4,3] + 1 * M[3,1] + 1 * M[1,1] = 0,
        3 * M[4,4] + 2 * M[5,4] + 1 * M[4,1] + 1 * M[2,1] = 0,
        ...

    In the standard form, we represent this condition by

        tr<C_(x,1,1), M> = 0,
        tr<C_(x,1,2), M> = 0,
        tr<C_(x,1,3), M> = 0,
        tr<C_(x,2,1), M> = 0,
        ...
        tr<C_(x,3,3), M> = 0,
        tr<C_(y,1,1), M> = 0,
        ...
        tr<C_(y,3,3), M> = 0,

        where C_(z, k, l) is the coefficient matrix corresponding to the condition that f_z * M'[k,l] = 0.

    So, in the standard form, the SDP problem for polynomial optimization with the first order condition is (we go back to the case that orderObj = 4)

        n = 1, p = 1, r_1 = 6,

        minimize        0 * x_1 + tr<C, X>

        by choosing     x_1 and X

        subject to      X[1,1] = 1,                                             The value of the zero-order monomial must be 1.
                        X[4,1] - X[2,2] = 0,                                    The monomial 'x2' appears multiple times in M, and their values must be equal.
                        X[5,1] - X[3,2] = 0,                                    The monomial 'xy' appears multiple times in M, and their values must be equal.
                        ...,
                        X[6,2] - X[5,3] = 0,                                    The monomial 'xy2' appears multiple times in M, and their values must be equal.
                        tr<C_(x,1,1), X> = 0,                                   gradient(f, x) * M[1,1] = 0,
                        tr<C_(y,1,1), X> = 0,                                   gradient(f, y) * M[1,1] = 0,
                        X is positive semidefinite.

    */

    /*

    Now we write the SDP formulation of the polynomial optimization problem
        in the structure that Mosek can understand.

    */

    /*

    We first count the number of constraints.

    For brevity of notation in the exposition, let
        dim  = dimMomentMatrix,
        dim' = dimLocalizingMatrix.

    First, there is a constraint that X[1,1] = 1.

    In addition, there are
        dim + choose(dim,2) = dim*(dim+1)/2     lower triangular entries in momentMatrix, and
        choose(nvar+order, nvar)                unique monomials.

    So there are 
        dim*(dim+1)/2 - choose(nvar+order, nvar)
    number of constraints that impose equivalence of monomial values.

    Also, for the first order conditions,
        we additionally have dim' + choose(dim',2) = dim'*(dim'+1)/2 constraints
        for each variable of the polynomial.

    Therefore, there are
        nvar * dim' * (dim'+1)/2
    number of constraints that impose first order conditions.

    To code the constraints in the structure that Mosek can understand, 
        we also need to compute the number of nonzero entries of the coefficient matrix
        appearing in the constraints.

    The number of nonzero entries are computed as follows:
        First, there is one nonzero entry in the constraint which is X[1,1] = 1.
        Second, there are two nonzero entries for each equivalence-of-monomial-values constraints.
        Third, there are _number_of_monomials(f_i)_ nonzero entries for each variable dimension and for each element in the localizing matrix.

    */

    // record nearest even degrees of momentMatrix and objective
    int evenOrderObj = 2 * floor( (orderObj(0) + 1) / 2 );
    int evenOrderMom = 2 * floor( (orderMom(0) + 1) / 2 );

    // compute the number of monomials up to _order_ -th order in _nvar_ variables
    int nMonomials = choose_cpp(static_cast<int>(nvar(0) + evenOrderMom), static_cast<int>(nvar(0)));

    // compute the dimension of momentMatrix. Note that it equals to the number of monomials up to as.integer( (_order_+1)/2 ).
    int dimMomentMatrix = choose_cpp(static_cast<int>(nvar(0) + evenOrderMom/2), static_cast<int>(nvar(0)));
    int dimLocalizingMatrix = choose_cpp(static_cast<int>(nvar(0) + evenOrderMom/2 - evenOrderObj/2), static_cast<int>(nvar(0)));

    // initialize containers for partial derivative list and its coefs and degrees
    List derivative;
    NumericVector derivCoefs;
    NumericMatrix derivDegrees;

    // compute lengths of gradient polynomials
    NumericVector lenGradient(static_cast<int>(nvar(0)));
    for(int d=0; d<nvar(0); d++) {
        derivative = gradientObj[d];
        derivCoefs = derivative["coefs"];
        lenGradient(d) = derivCoefs.size();
    }

    // number of constraints
    int entriesInMomentMatrix = dimMomentMatrix * (dimMomentMatrix + 1) / 2;
    int entriesInLocalizingMatrix = dimLocalizingMatrix * (dimLocalizingMatrix + 1) / 2;
    int nConstraints = 1 + 
                       entriesInMomentMatrix - nMonomials + 
                       nvar(0) * entriesInLocalizingMatrix;

    // number of nonzero coefficient matrix entries in the constraints
    int nNonzeroEntries = 1 + 
                          2 * ( entriesInMomentMatrix - nMonomials ) +
                          sum(lenGradient * entriesInLocalizingMatrix);

    // dimension of the scalar choice variable (x_1 in the standard form) is 1, so we don't define it.
    ///// int nScalarVariable = 1; /////

    /*

    Now we construct Mosek's SDP program structure.

    We first define each piece of the structure and then combine them as a list at the end.

    */

    // optimization sense, either "max" or "min". Left blank in this function
    CharacterVector sense("assign max or min, depending on the problem.");
    // coefficient on the scalar variable in the objective, which is set to zero by default
    NumericVector c(1);
    // coefficient on the scalar variable in the constraints, which are set to zeros by default
    NumericMatrix A(nConstraints,1);
    // lower and upper bounds for the constraints. The entries in the first column are one (which represents X[1,1] = 1), otherwise zero.
    NumericMatrix bc(2,nConstraints);
    bc(0,0) = 1; bc(1,0) = 1;
    // lower and upper bounds for the scalar variable, which are zeros
    NumericMatrix bx(2,1);

    // dimension of the matrix variable
    NumericVector bardim = NumericVector::create(dimMomentMatrix);

    // coefficients on the matrix variable in the objective.
    // this part is left blank since it is problem-specific and not a "skeleton".
    CharacterVector barc("fill in the coefficients in the sparseMatrix triplet form.");

    /*

    coefficients on the matrix variable in the constraints, barA.

    */

    // we define each entry of barA separately.
    NumericVector barA_i(nNonzeroEntries);
    NumericVector barA_j(nNonzeroEntries);
    NumericVector barA_k(nNonzeroEntries);
    NumericVector barA_l(nNonzeroEntries);
    NumericVector barA_v(nNonzeroEntries);

    // first write the X[1,1] = 1 constraint
    barA_i(0) = 1;
    barA_j(0) = 1;
    barA_k(0) = 1;
    barA_l(0) = 1;
    barA_v(0) = 1;

    // define iterators and counters needed to enumerate the rest of the constraints
    int barAInd = 1;   // position in the bara vectors
    int iInd = 2;      // counter for the constraint number
    int j, k, l;
    NumericVector momentId(1);
    NumericVector momentIdPosition(1);

    /*

    write equivalence of the monomial values

    */

    // enumerate unique elements in momentMatrix
    NumericVector momentVector = momentMatrixSparse(_,2);
    NumericVector uniqueMomentVector = unique(momentVector);

    // impose equivalence for each entry in uniqueMomentVector
    for(k=0; k<uniqueMomentVector.size(); k++) {

        // read moment ID
        momentId = uniqueMomentVector(k);
        // record the earliest position of momentID in momentVector
        momentIdPosition = match(momentId, momentVector) - 1;

        // if moment appears again in the remaining sequence of the momentVector, record the equivalence constraint
        // search for additional appearance if the first appearance is not the last entry of momentVector
        if(momentIdPosition(0) < momentVector.size()) {
            for(l=momentIdPosition(0)+1; l<momentVector.size(); l++) {
                // if there is an additional appearance, record the constraint
                if(momentId(0) == momentVector(l)) {
                    // code the constraint
                    barA_i(barAInd) = iInd;
                    barA_j(barAInd) = 1;
                    barA_k(barAInd) = momentMatrixSparse(momentIdPosition(0),0);
                    barA_l(barAInd) = momentMatrixSparse(momentIdPosition(0),1);
                    // if off-diagonal, deflate the coefficient by half since the same is also copied to the upper-triangular part in Mosek
                    if(barA_k(barAInd) == barA_l(barAInd)) {
                        barA_v(barAInd) = 1.0;
                    } else {
                        barA_v(barAInd) = 0.5;
                    }
                    
                    barA_i(barAInd+1) = iInd;
                    barA_j(barAInd+1) = 1;
                    barA_k(barAInd+1) = momentMatrixSparse(l,0);
                    barA_l(barAInd+1) = momentMatrixSparse(l,1);
                    // if off-diagonal, deflate the coefficient by half since the same is also copied to the upper-triangular part in Mosek
                    if(barA_k(barAInd+1) == barA_l(barAInd+1)) {
                        barA_v(barAInd+1) = -1.0;
                    } else {
                        barA_v(barAInd+1) = -0.5;
                    }
                    iInd++;
                    barAInd = barAInd + 2;
                }
            }
        }
    }

    /*

    write first order conditions

    */

    // write constraint for each variable dimension
    for(j=0; j<nvar(0); j++) {
        // read gradient information
        derivative   = gradientObj[j];
        derivCoefs   = derivative["coefs"];
        derivDegrees = SEXP(derivative["degrees"]);
        // write constraint for each entry in the upper triangular part of the momentMatrix
        for(k=0; k<momentMatrixSparse.nrow(); k++) {
            // if momentMatrix entry belongs to localizing matrix
            if(momentMatrixSparse(k,0) <= dimLocalizingMatrix) {
                // record for each monomial in the gradient polynomial
                for(l=0; l<lenGradient(j); l++) {
                    // get monomial ID for the l-th monomial of partial derivative j multiplied with k-th sparse-entry of localizing matrix
                    momentId = degreeToMonomial_cpp(derivDegrees(l,_), monomialPrimes) * momentVector(k);
                    // find its position in momentMatrixSparse
                    momentIdPosition = match(momentId, momentVector) - 1;
                    // record the coefficient
                    barA_i(barAInd) = iInd;
                    barA_j(barAInd) = 1;
                    barA_k(barAInd) = momentMatrixSparse(momentIdPosition(0),0);
                    barA_l(barAInd) = momentMatrixSparse(momentIdPosition(0),1);
                    // if off-diagonal, deflate the coefficient by half since the same is also copied to the upper-triangular part in Mosek
                    if(barA_k(barAInd) == barA_l(barAInd)) { barA_v(barAInd) = derivCoefs(l); } else { barA_v(barAInd) = 0.5 * derivCoefs(l); }
                    barAInd++;
                }
                iInd++;
            }
        }
    }

    // return the skeleton
    return List::create(_["sense"] = sense,
                        _["c"] = c,
                        _["A"] = A,
                        _["bc"] = bc,
                        _["bx"] = bx,
                        _["bardim"] = bardim,
                        _["barc"] = barc,
                        _["barA"] = List::create(_["i"] = barA_i,
                                                 _["j"] = barA_j,
                                                 _["k"] = barA_k,
                                                 _["l"] = barA_l,
                                                 _["v"] = barA_v));

}