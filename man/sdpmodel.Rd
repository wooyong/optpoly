\name{sdpmodel}
\alias{sdpmodel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create SDP Model Object for Global Optimization
}
\description{
Create a semidefinite programming model for global optimization of a polynomial
}
\usage{
sdpmodel(sense, coefs, degrees, opt = NULL)
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
This function creates an R object that describes the semidefinite programming model proposed by Lasserre (2001) and Nie et al (2006) for global optimization of a polynomial. The model is described according to Rmosek's modeling rule, and it can directly be supplied to Rmosek's \code{mosek} function.
}
\value{
A hierarchy of semidefinite programs where each step contains:
\item{model}{a list specifying semidefinite programming model for global optimization of a polynomial}
\item{monomialSystem}{an output of \code{createMonomialVector} used to describe \code{model}}
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
# 4x1^2 - 2.1x1^4 + (1/3)x1^6 + x1x2 - 4x2^2 + 4x2^4
cf=c(4,-2.1,1/3,1,-4,4)
dg=matrix(c(2,0,4,0,6,0,1,1,0,2,0,4), ncol=2, byrow=TRUE)

sdpmodel("min", cf, dg)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
