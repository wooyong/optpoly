\name{createMomentMatrixSparse}
\alias{createMomentMatrixSparse}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create Moment Matrix
}
\description{
Create symbolic moment matrix of the semidefinite program in sparse triplet form
}
\usage{
createMomentMatrixSparse(monomialVector)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{monomialVector}{
%%     ~~Describe \code{monomialVector} here~~
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
function (monomialVector)
{
    momentMatrix = createMomentMatrix(monomialVector)
    momentMatrixSparse = data.frame(k = as.vector(row(momentMatrix)),
        l = as.vector(col(momentMatrix)), v = as.vector(momentMatrix))
    momentMatrixSparse = subset(momentMatrixSparse, k >= l)
    return(momentMatrixSparse)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
