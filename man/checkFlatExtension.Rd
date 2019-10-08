\name{checkFlatExtension}
\alias{checkFlatExtension}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Check Flat Extension Condition
}
\description{
Check flat extension rank condition of the semidefinite solution
}
\usage{
checkFlatExtension(optimalMomentMatrix, varDim, largeOrder, smallOrder, ereltol)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{optimalMomentMatrix}{
%%     ~~Describe \code{optimalMomentMatrix} here~~
}
  \item{varDim}{
%%     ~~Describe \code{varDim} here~~
}
  \item{largeOrder}{
%%     ~~Describe \code{largeOrder} here~~
}
  \item{smallOrder}{
%%     ~~Describe \code{smallOrder} here~~
}
  \item{ereltol}{
%%     ~~Describe \code{ereltol} here~~
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
function (optimalMomentMatrix, varDim, largeOrder, smallOrder,
    ereltol)
{
    largeMonomialOrder = as.integer((largeOrder + 1)/2)
    smallMonomialOrder = as.integer((smallOrder + 1)/2)
    largeMonomialSystem = createMonomialVector(varDim, largeMonomialOrder)
    smallMonomialSystem = createMonomialVector(varDim, smallMonomialOrder)
    largeMonomialLen = length(largeMonomialSystem$vec)
    smallMonomialLen = length(smallMonomialSystem$vec)
    largeMatrix = as.matrix(optimalMomentMatrix[1:largeMonomialLen,
        1:largeMonomialLen])
    smallMatrix = as.matrix(optimalMomentMatrix[1:smallMonomialLen,
        1:smallMonomialLen])
    largeRank = rankSymmetricMatrix(largeMatrix, ereltol)
    smallRank = rankSymmetricMatrix(smallMatrix, ereltol)
    return(list(largeRank = largeRank, smallRank = smallRank,
        certificate = (largeRank == smallRank)))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
