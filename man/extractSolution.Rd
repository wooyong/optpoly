\name{extractSolution}
\alias{extractSolution}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extract Solution
}
\description{
Extract global optimizers of a polynomial from the semidefinite solution
}
\usage{
extractSolution(sol, points = NULL, vreltol = 1e-10, ereltol = 0.001)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sol}{
output of the function \code{optpoly}
}
  \item{points}{
number of optimizers, set to \code{sol$smallRank} if missing
}
  \item{vreltol}{
tolerance for regarding a value smaller than \code{vreltol} in absolute value as zero
}
  \item{ereltol}{
tolerance for regarding a relative drop in the eigenvalue sequence as a drop to zero
}
}
\details{
This function extracts optimal solution of the polynomial optimization problem given the output of \code{optpoly}. In theory, the number of solutions equals to the rank of the optimal moment matrix given by the semidefinite program. However, since it is hard to numerically check the rank of the matrix, the user may also directly supply the number of optimizers to \code{points}. Sometimes The optimal solution can be found by trying different values to \code{points}.
}
\value{
A matrix where each row represents a solution to the polynomial optimization problem.
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

sol=optpoly("min", cf, dg, list(radius=3))
extractSolution(sol, points=2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
