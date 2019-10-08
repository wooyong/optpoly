\name{createMonomialVector}
\alias{createMonomialVector}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create Monomial Vector
}
\description{
Create symbolic monomial vector used for describing semidefinite programs
}
\usage{
createMonomialVector(varDim, order)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{varDim}{
%%     ~~Describe \code{varDim} here~~
}
  \item{order}{
%%     ~~Describe \code{order} here~~
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
function (varDim, order)
{
    lengthMatrix = matrix(0, nrow = max(order, 1), ncol = varDim)
    lengthMatrix[1, ] = rep(1, varDim)
    if (order >= 2) {
        for (r in 2:order) {
            for (v in 1:varDim) {
                lengthMatrix[r, v] = sum(lengthMatrix[r - 1,
                  v:varDim])
            }
        }
    }
    primes = getPrimes(as.integer(2 + 2 * varDim * log(varDim)))
    monomialVectorList = vector("list", 1 + order)
    monomialVectorList[[1]] = 1
    if (order >= 1) {
        monomialVectorList[[2]] = primes[1:varDim]
    }
    if (order >= 2) {
        for (r in 2:order) {
            for (v in 1:varDim) {
                if (v == 1) {
                  idVec = primes[v] * monomialVectorList[[r]]
                }
                else {
                  idVec = primes[v] * monomialVectorList[[r]][(1 +
                    sum(lengthMatrix[r - 1, 1:(v - 1)])):sum(lengthMatrix[r -
                    1, ])]
                }
                monomialVectorList[[1 + r]] = c(monomialVectorList[[1 +
                  r]], idVec)
            }
        }
    }
    monomialVector = unlist(monomialVectorList, use.names = FALSE)
    return(list(vec = monomialVector, primes = primes[1:varDim]))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
