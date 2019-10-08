\name{optpoly}
\alias{optpoly}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Global Optimization of a Polynomial
}
\description{
Globally optimize a polynomial using semidefinite relaxation approaches proposed by Lasserre (2001) and Nie et al (2006).
}
\usage{
optpoly(sense, coefs, degrees, opt = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sense}{
whether to minimize or maximize the polynomial, either \code{"min"} or \code{"max"}.
}
  \item{coefs}{
a vector of coefficients, \code{length(coefs)} must equal to \code{nrow(degrees)}.
}
  \item{degrees}{
a matrix of polynomial degrees, each row represents powers of the monomial for each entry of \code{coefs}.
}
  \item{opt}{
a list of options.
}
}
\details{
This function solves global polynomial optimization with or without the constraint that the domain of the polynomial is in the sphere of the specified radius in \code{opt}. When unconstrained, this function implements Nie et al (2006)'s optimization method based on the gradient ideal. It solves a hierarchy of semidefinite programs (or a sequence of semidefinite programs) which converges to the exact solution. When constrained, this function implements Lasserre (2001)'s optimization method. It solves a hierarchy of semidefinite programs (or a sequence of semidefinite programs) which converges to the exact solution. The user can specify the following options:
\describe{
\item{\code{radius}}{a positive real number. If supplied, the domain of the polynomial is restricted to be in the sphere of length \code{radius}. If omitted, the domain is assumed to be the entire Euclidean space.}
\item{\code{hierarchy}}{an integer that represents the length of the hierarchy of the semidefinite relaxations. Default is \code{1}.}
\item{\code{multithread}}{an indicator for whether to allow for multithreading and utilize multiple CPU cores. Default is \code{FALSE}.}
\item{\code{verbose}}{an integer variable that is passed to the SDP solver Mosek. If positive, Mosek will print its output on the R console. Default is \code{0}.}
\item{\code{ereltol}}{a positive real number used in checking the flat extension condition for certificate of the exact solution. It specifies how much relative drop in the eigenvalue sequence is interpreted as a drop to zero. Default is \code{1e-03}.}
}
}
\value{
\item{objective_primal}{primal objective of the semidefinite program}
\item{objective_dual}{dual objective of the semidefinite program}
\item{sdpstatus}{status of the semidefinite program}
\item{solstatus}{status of the solution}
\item{momentmatrix}{optimal moment matrix}
\item{varDim}{dimension of the polynomial interpreted from the input}
\item{order}{order of the polynomial interpreted from the input}
\item{largeRank}{rank of the "large" matrix in the flat extension condition}
\item{smallRank}{rank of the "small" matrix in the flat extension condition}
\item{certificate}{certificate for exact solution, \code{TRUE} if \code{largeRank == smallRank}}
\item{hierarchy}{number of hierarchies solved to obtain current output}
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
require(Rmosek)

# 4x1^2 - 2.1x1^4 + (1/3)x1^6 + x1x2 - 4x2^2 + 4x2^4
cf=c(4,-2.1,1/3,1,-4,4)
dg=matrix(c(2,0,4,0,6,0,1,1,0,2,0,4), ncol=2, byrow=TRUE)

optpoly("min", cf, dg, list(radius=3))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
