\name{optpoly}
\alias{optpoly}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Global Optimization of a Polynomial
}
\description{
Globally optimize a polynomial using semidefinite relaxation approach proposed by Lasserre (2001) and Nie et al (2006).
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
a matrix of polynomial degrees, each row represents powers of monomial for each entry of \code{coefs}.
}
  \item{opt}{
a list of options.
}
}
\details{
This function solves global polynomial optimization with or without the constraint that the domain of the polynomial is a rectangle or a sphere, specified in \code{opt}. When unconstrained, this function implements Nie et al (2006)'s optimization method based on gradient ideals. When constrained, this function implements Lasserre (2001)'s optimization method. In both cases, a hierarchy of semidefinite programs (SDP), or a sequence of SDPs, is solved, which converges to exact solution. The user can specify the following options:
\describe{
\item{\code{bounds}}{a matrix of \code{2} rows and \code{ncol(degrees)} columns. If supplied, the domain of the polynomial is restricted to be a rectangle with the first row being the lower bound and the second row being the upper bound. If omitted, the domain is assumed to be the entire Euclidean space.}
\item{\code{radius}}{a positive real number. If supplied, the domain of the polynomial is restricted to be the sphere of length \code{radius} centered at the origin. If omitted, the domain is assumed to be the entire Euclidean space.}
\item{\code{hierarchy}}{an integer. The length of hierarchy of semidefinite relaxations. Default is \code{1}.}
\item{\code{multithread}}{an indicator for whether to allow for multithreading to solve SDPs. Default is \code{FALSE}.}
\item{\code{verbose}}{an integer that is passed to the SDP solver Mosek. If positive, Mosek will print its output on the R console. Default is \code{0}.}
\item{\code{fabstol}}{absolute tolerance for function value used for certificate of optimality. Default is \code{1e-10}.}
\item{\code{freltol}}{relative tolerance for function value used for certificate of optimality. Default is \code{1e-06}.}
\item{\code{vabstol}}{absolute tolerance for regarding a value smaller than \code{vabstol} in absolute value as zero. Default is \code{1e-10}.}
\item{\code{ereltol}}{relative tolerance for regarding a drop in the eigenvalue sequence as a drop to zero. It is used for checking the rank condition for certificate of optimality. Default is \code{1e-03}.}
}
}
\value{
\item{objective_primal}{primal objective value of the SDP}
\item{objective_dual}{dual objective value of the SDP}
\item{sdpstatus}{status of the SDP}
\item{solstatus}{status of the SDP solution}
\item{momentmatrix}{optimal moment matrix}
\item{varDim}{dimension of the polynomial interpreted from the input}
\item{order}{order of the polynomial interpreted from the input}
\item{certificate}{certificate of optimality for exact solution}
\item{rank}{number of optimal solutions, if \code{certificate = TRUE}}
\item{hierarchy}{number of iterations performed within hierarchy to obtain current output}
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

# minimization over R^2 space
optpoly("min", cf, dg)

# minimization over a circle of radius 3 centered at the origin
optpoly("min", cf, dg, list(radius=3))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
