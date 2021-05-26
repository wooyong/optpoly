\name{solvepoly}
\alias{solvepoly}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Solve a System of Polynomial Equations
}
\description{
Solve a system of polynomial equations using semidefinite relaxation approach proposed by Lasserre (2001) and Nie et al (2006).
}
\usage{
solvepoly(equationsList, opt = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{equationsList}{
a list of polynomial equations.
}
  \item{opt}{
a list of options.
}
}
\details{
This function solves a system of polynomial equations using an algorithm that solves hierarchy of semidefinite programs (SDP). This algorithm only searches for real-valued solutions, which is computationally cheaper than algorithms that search for all complex solutions. The system can be solved over (i) entire Euclidean space (unconstrained), (ii) a sphere centered at origin or (iii) a rectangle. This function uses unconstrained domain by default, and the user can specify sphere or rectangular domain in \code{opt} to use constrained domain.

A system of polynomial equations is specified as an unlabeled list in \code{equationsList}. Each element of a list represents one polynomial equation, where the right-hand side is zero and the left-hand side is a polynomial represented by:
\describe{
\item{coefs}{
a vector of polynomial coefficients.
}
\item{degrees}{
a matrix of polynomial degrees, each row representing powers of monomial for each entry of \code{coefs}.
}
}

To solve a system of polynomial equations, the algorithm solves an optimization problem where it minimizes a zero function subject to the constraint that the system of polynomial equations hold. Similar to the function \code{optpoly}, this amounts to solving a hierarchy of the SDPs, or a sequence of the SDPs, which converges to an exact solution. If the optimization gives an exact solution, the user can take the output to the function \code{extractSolution} and obtain solutions to the system of equations.

This function uses modified version of \code{optpoly} for the optimization, and so it accepts the same set of options in \code{opt}. The user can specify the following options:
\describe{
\item{\code{bounds}}{a matrix of \code{2} rows and \code{ncol(degrees)} columns. If supplied, optimization is performed over a rectangle with the lower bound being the first row and the upper bound being the second row. If both \code{bounds} and \code{radius} are omitted, optimization is performed over the entire Euclidean space.}
\item{\code{radius}}{a positive real number. If supplied, optimization is performed over a sphere of length \code{radius} centered at the origin. If both \code{bounds} and \code{radius} are omitted, optimization is performed over the entire Euclidean space.}
\item{\code{hierarchy}}{an integer. The length of hierarchy of semidefinite relaxations. Default is \code{1}.}
\item{\code{multithread}}{an indicator for using multithreading for solving SDPs. Default is \code{FALSE}.}
\item{\code{verbose}}{an integer. If positive, Mosek will print its output on the R console (see Rmosek documentation). Default is \code{0}.}
\item{\code{fabstol}}{absolute tolerance for polynomial value used for certificate of optimality. Default is \code{1e-10}.}
\item{\code{freltol}}{relative tolerance for polynomial value used for certificate of optimality. Default is \code{1e-06}.}
\item{\code{vabstol}}{absolute tolerance for regarding a value smaller than \code{vabstol} in absolute value as zero. Default is \code{1e-10}.}
\item{\code{ereltol}}{relative tolerance for regarding a drop in the eigenvalue sequence as a drop to zero. It is used for checking the rank condition for certificate of optimality. Default is \code{1e-03}.}
}
}
\value{
\item{objective_primal}{primal objective value of the SDP}
\item{objective_dual}{dual objective value of the SDP}
\item{sdpstatus}{status of the SDP, passed from Rmosek}
\item{solstatus}{status of the SDP solution, passed from Rmosek}
\item{momentmatrix}{optimal moment matrix}
\item{varDim}{dimension of the supplied polynomial}
\item{order}{degree of the supplied polynomial}
\item{certificate}{certificate of optimality}
\item{rank}{number of optimal solutions, valid only if \code{certificate = TRUE}}
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

# define system of equations
eqList = list(
    # x^2 + y^2 - 20 = 0
    list(degrees = matrix(c(2,0,  # x^2
                            0,2,  # y^2
                            0,0), # constant
                            nrow=3, ncol=2, byrow=TRUE),
         coefs   = c(1, 1, -20)),
    # x + y - 6 = 0
    list(degrees = matrix(c(1,0,
                            0,1,
                            0,0),
                            nrow=3, ncol=2, byrow=TRUE),
        coefs    = c(1, 1, -6))
)

# solve system of equations via optimization
sol = solvepoly(eqList, opt=list(hierarchy=2))

# extract solution from exact optimization result
extractSolution(sol)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
