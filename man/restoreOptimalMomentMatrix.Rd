\name{restoreOptimalMomentMatrix}
\alias{restoreOptimalMomentMatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Restore Optimal Moment Matrix
}
\description{
Restore optimal moment matrix from the semidefinite solution
}
\usage{
restoreOptimalMomentMatrix(monomialSystem, optimalMomentVector)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{monomialSystem}{
%%     ~~Describe \code{monomialSystem} here~~
}
  \item{optimalMomentVector}{
%%     ~~Describe \code{optimalMomentVector} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (monomialSystem, optimalMomentVector)
{
    dimMomentMatrix = length(monomialSystem$vec)
    optimalMomentMatrix = matrix(0, nrow = dimMomentMatrix, ncol = dimMomentMatrix)
    optimalMomentMatrix[row(optimalMomentMatrix) >= col(optimalMomentMatrix)] = optimalMomentVector
    optimalMomentMatrix = optimalMomentMatrix + t(optimalMomentMatrix)
    diag(optimalMomentMatrix) = diag(optimalMomentMatrix)/2
    return(optimalMomentMatrix)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
