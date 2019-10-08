\name{optmanypoly}
\alias{optmanypoly}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Global Optimization of Many Polynomials
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
optmanypoly(sense, coefMatrix, varDim, orderObj, orderMom, constrained, radius = NULL, multithread, verbose)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sense}{
%%     ~~Describe \code{sense} here~~
}
  \item{coefMatrix}{
%%     ~~Describe \code{coefMatrix} here~~
}
  \item{varDim}{
%%     ~~Describe \code{varDim} here~~
}
  \item{orderObj}{
%%     ~~Describe \code{orderObj} here~~
}
  \item{orderMom}{
%%     ~~Describe \code{orderMom} here~~
}
  \item{constrained}{
%%     ~~Describe \code{constrained} here~~
}
  \item{radius}{
%%     ~~Describe \code{radius} here~~
}
  \item{multithread}{
%%     ~~Describe \code{multithread} here~~
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
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
function (sense, coefMatrix, varDim, orderObj, orderMom, constrained,
    radius = NULL, multithread, verbose)
{
    N = nrow(coefMatrix)
    monomialSystemObj = createMonomialVector(varDim, as.integer((orderObj +
        1)/2))
    momentMatrixSparseObj = as.matrix(createMomentMatrixSparse(monomialSystemObj$vec))
    monomialSystemMom = createMonomialVector(varDim, as.integer((orderMom +
        1)/2))
    momentMatrixSparseMom = as.matrix(createMomentMatrixSparse(monomialSystemMom$vec))
    degrees = monomialsToDegrees_cpp(momentMatrixSparseObj[,
        3], monomialSystemObj$primes)
    if (constrained == 0) {
        res = foreach(i = 1:N, .combine = "rbind") \%dopar\% {
            grad = computeGradient(coefMatrix[i, ], degrees,
                monomialSystemObj$primes)
            model = createMosekSdpModelSkeletonWithGradientIdeals_cpp(varDim,
                orderObj, orderMom, grad, monomialSystemMom$primes,
                momentMatrixSparseMom)
            model = addLabelsToSkeleton(model, sense, multithread)
            model$barc = createMosekSdpCoefficientMatrixFromMonomials_cpp(coefMatrix[i,
                ], momentMatrixSparseObj[, 3], momentMatrixSparseObj)
            mosekSol = mosek(model, opts = list(verbose = verbose,
                soldetail = 1))
            data.frame(objective_primal = mosekSol$sol$itr$pobjval,
                objective_primal = mosekSol$sol$itr$pobjval,
                sdpstatus = mosekSol$response$msg, solstatus = mosekSol$sol$itr$solsta,
                matrix(mosekSol$sol$itr$barx[[1]], nrow = 1,
                  ncol = nrow(momentMatrixSparseMom)))
        }
    }
    else {
        model = createMosekSdpModelSkeleton_cpp(varDim, orderMom,
            constrained, radius, monomialSystemMom$primes, momentMatrixSparseMom)
        model = addLabelsToSkeleton(model, sense, multithread)
        res = foreach(i = 1:N, .combine = "rbind") \%dopar\% {
            model$barc = createMosekSdpCoefficientMatrixFromMonomials_cpp(coefMatrix[i,
                ], momentMatrixSparseObj[, 3], momentMatrixSparseObj)
            mosekSol = mosek(model, opts = list(verbose = verbose,
                soldetail = 1))
            data.frame(objective_primal = mosekSol$sol$itr$pobjval,
                objective_primal = mosekSol$sol$itr$pobjval,
                sdpstatus = mosekSol$response$msg, solstatus = mosekSol$sol$itr$solsta,
                matrix(mosekSol$sol$itr$barx[[1]], nrow = 1,
                  ncol = nrow(momentMatrixSparseMom)))
        }
    }
    return(list(objective_primal = res$objective_primal, objective_dual = res$objective_dual,
        sdpstatus = factor(res$sdpstatus), solstatus = factor(res$solstatus),
        varDim = varDim, orderobj = orderObj, orderMom = orderMom,
        moment_vectors = res[, -c(1:4)]))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
