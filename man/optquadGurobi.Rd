\name{optquadGurobi}
\alias{optquadGurobi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Globally Optimize a Quadratic Polynomial
}
\description{
Globally optimize a quadratic polynomial using quadratic programming (QP) solver in Gurobi
}
\usage{
optquadGurobi(sense, coefs, degrees, constraints = NULL, opt = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sense}{
to minimize or maximize the polynomial, either \code{"min"} or \code{"max"}.
}
  \item{coefs}{
a vector of polynomial coefficients.
}
  \item{degrees}{
a matrix of polynomial degrees, each row representing powers of monomial for each entry of \code{coefs}.
}
  \item{constraints}{
a list of constraint configurations.
}
  \item{opt}{
a list of options.
}
}
\details{
This function optimizes a quadratic polynomial using QP solver in Gurobi. This function solves unconstrained optimization by default, and the user can specify constraints to perform constrained optimization. This function accepts constraints of the form \code{Ax (<=, =, >=) rhs} and \code{lb <= x <= ub} where \code{A} is a matrix, \code{x} is a choice variable vector and \code{rhs,lb,ub} are lower and upper bound vectors. The user can specify them in \code{constraints} by a list containing (some can be omitted if there is no constraint of the corresponding form):
\describe{
\item{\code{A}}{a coefficient matrix.}
\item{\code{rhs}}{a right-hand side vector.}
\item{\code{sense}}{a vector of strings, which are either \code{"<="}, \code{"="} or \code{">="}.}
\item{\code{lb}}{a vector of lower bounds.}
\item{\code{ub}}{a vector of upper bounds.}
}
The user can also specify the following options in \code{opt}:
\describe{
\item{\code{Threads}}{the number of threads used for solving QP. Default is \code{1}.}
\item{\code{OutputFlag}}{an integer. If \code{1}, Gurobi will print its output. If \code{0}, it will not (see Gurobi documentation). Default is \code{1}.}
}
}
\value{
\item{objective}{objective value of the QP}
\item{qpstatus}{status of the QP, passed from Gurobi}
\item{solution}{optimal solution}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Wooyong Lee
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
require(gurobi)

# x^2 + y^2 + 2xy + 2x + 2y + 2
cf=c(1,1,2,2,2,2)
dg=matrix(c(2,0,0,2,1,1,1,0,0,1,0,0), nrow=6, ncol=2, byrow=TRUE)

# minimization over x+y <= 1 and -1 <= x,y <=1
A=matrix(c(1,1), nrow=1, ncol=2)
rhs=1
sense="<="
lb=c(-1,1)
ub=c(-1,1)
optquadGurobi("min", cf, dg, constraints=list(A=A, rhs=rhs, sense=sense, lb=lb, ub=ub))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
