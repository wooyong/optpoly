\name{optquad}
\alias{optquad}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Globally Optimize a Quadratic Polynomial
}
\description{
Globally optimize a quadratic polynomial using quadratic programming (QP) solver in Mosek
}
\usage{
optquad(sense, coefs, degrees, constraints = NULL, opt = NULL)
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
This function optimizes a quadratic polynomial using QP solver in Mosek. This function solves unconstrained optimization by default, and the user can specify constraints to perform constrained optimization. This function accepts constraints of the form \code{blc <= Ax <= buc} and \code{blx <= x <= bux} where \code{A} is a matrix, \code{x} is a choice variable vector and \code{blc,buc,blx,bux} are lower and upper bound vectors. The user can specify them in \code{constraints} by a list containing (some can be omitted if there is no constraint of the corresponding form):
\describe{
\item{\code{A}}{a coefficient matrix.}
\item{\code{bc}}{a matrix of two rows, the first row being \code{blc} and the second row being \code{buc}.}
\item{\code{bx}}{a matrix of two rows, the first row being \code{blx} and the second row being \code{bux}.}
}
The user can also specify the following options in \code{opt}:
\describe{
\item{\code{multithread}}{an indicator for using multithreading for solving QPs. Default is \code{FALSE}.}
\item{\code{verbose}}{an integer. If positive, Mosek will print its output on the R console (see Rmosek documentation). Default is \code{0}.}
}
}
\value{
\item{objective_primal}{primal objective value of the QP}
\item{objective_dual}{dual objective value of the QP}
\item{qpstatus}{status of the QP, passed from Rmosek}
\item{solstatus}{status of the QP solution, passed from Rmosek}
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
require(Rmosek)

# x^2 + y^2 + 2xy + 2x + 2y + 2
cf=c(1,1,2,2,2,2)
dg=matrix(c(2,0,0,2,1,1,1,0,0,1,0,0), nrow=6, ncol=2, byrow=TRUE)

# minimization over x+y <= 1 and -1 <= x,y <=1
A=matrix(c(1,1), nrow=1, ncol=2)
bc=matrix(c(-Inf,1), nrow=2, ncol=1)
bx=matrix(c(-1,-1,1,1), nrow=2, ncol=2, byrow=TRUE)
optquad("min", cf, dg, constraints=list(A=A, bc=bc, bx=bx))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
