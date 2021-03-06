\name{createMonomialVector}
\alias{createMonomialVector}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create Monomial Vector
}
\description{
Create symbolic vector of monomials.
}
\usage{
createMonomialVector(varDim, order)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{varDim}{
dimension of variables.
}
  \item{order}{
maximum degree of monomials.
}
}
\details{
This function creates a complete list of monomials in \code{varDim} variables of degree at most \code{order}. To represent monomials numerically, I map variables into prime numbers. For example, to represent a monomial in two variables, I map (x1, x2) -> (2, 3). Then a monomial x1^2*x2 is represented by 2^2*3 = 12.
}
\value{
\item{vec}{a vector containing all monomials in \code{varDim} variables of degree at most \code{order}}
\item{primes}{a vector of prime numbers representing each variable}
\item{order}{maximum degree of monomials}
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
# list of monomials in two dimensions of degree at most 2
createMonomialVector(2, 2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
