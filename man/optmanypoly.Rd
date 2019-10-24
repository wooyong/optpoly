\name{optmanypoly}
\alias{optmanypoly}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Global Optimization of Many Polynomials
}
\description{
Globally optimize many polynomials of same dimension and same degree.
}
\usage{
optmanypoly(sense, coefMatrix, varDim, orderObj, orderMom, constrained, radius = NULL, multithread, verbose)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sense}{
whether to minimize or maximize the polynomial, either \code{"min"} or \code{"max"}.
}
  \item{coefMatrix}{
a matrix of coefficients, each row representing a polynomial.
}
  \item{varDim}{
dimension of variables.
}
  \item{orderObj}{
degree of polynomials in \code{coefMatrix}.
}
  \item{orderMom}{
degree of semidefinite relaxation problems.
}
  \item{constrained}{
\code{0} if optimization is over R^n space, \code{1} if optimization is over a sphere of length \code{radius}.
}
  \item{radius}{
a positive real number. If supplied and \code{constrained = 1}, the domain of the polynomial is restricted to be in the sphere of length \code{radius}.
}
  \item{multithread}{
an indicator for whether to allow for multithreading to solve SDPs.
}
  \item{verbose}{
an integer that is passed to the SDP solver Mosek. If positive, Mosek will print its output on the R console.
}
}
\details{
This function solves global optimization problems of many polynomials of same dimension and degree, using the same level of semidefinite programming (SDP) relaxation. Each row of \code{coefMatrix} represents a polynomial, and the coefficients must be arranged in the same order as \code{createMomentMatrixSparse(varDim, orderMom)}. This function does not check certificate of optimality.
}
\value{
\item{objective_primal}{a vector of primal objectives of the SDPs}
\item{objective_dual}{a vector of dual objective of the SDPs}
\item{sdpstatus}{a vector of status of the SDPs}
\item{solstatus}{a vector of status of the SDP solutions}
\item{varDim}{dimension of polynomials}
\item{orderObj}{degree of polynomials in \code{coefMatrix}}
\item{orderMom}{degree of semidefinite relaxation problems}
\item{moment_vectors}{a matrix in which each row represents a moment matrix arranged in the same order as \code{createMomentMatrixSparse(varDim, orderMom)}}
}
\references{
Lasserre, Jean B (2001). Global optimization with polynomials and the problem of moments, \emph{SIAM Journal on optimization}, \strong{11}, 796-817.

Nie et al (2006). Minimizing polynomials via sum of squares over the gradient ideal, \emph{Mathematical programming}, \strong{106}, 587-606.

Lasserre, Jean B (2015). \emph{An introduction to polynomial and semi-algebraic optimization}, Cambridge University Press.
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

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
