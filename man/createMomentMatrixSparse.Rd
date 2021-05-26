\name{createMomentMatrixSparse}
\alias{createMomentMatrixSparse}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create Moment Matrix
}
\description{
Create symbolic moment matrix in sparse triplet form.
}
\usage{
createMomentMatrixSparse(monomialVector)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{monomialVector}{
a vector of monomials produced by \code{createMonomialVector$vec}.
}
}
\details{
This function creates a moment matrix, namely the outer product of \code{monomialvector}, in sparse triplet form. The maximum degree of monomials appearing in the moment matrix is twice of the maximum degree of \code{monomialvector}.
}
\value{
a moment matrix expressed in sparse triplet form.
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
monomialVector = createMonomialVector(2, 2)
# create outer product of monomialvector in sparse triplet form
createMomentMatrixSparse(monomialVector$vec)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
