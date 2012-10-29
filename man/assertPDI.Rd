\name{assertPDI}
\alias{assertPDI}
\title{Assert that a Matrix is an Element of P(D, I)}
\description{\code{assertPDI} checks if the matrix \code{M} is an element of \code{P(D, I)}. If not, the function produces an error.}
\usage{
assertPDI(M, D, I)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{M}{The matrix for which to check if it is an element of \code{P(D, I)}.}
  \item{D}{The adjacency matrix of the normal linear ADG model's independence graph sorted according to its topological order.}
  \item{I}{The index list \code{I} as returned from the \code{nlADG} class.}
}
\references{Andersson, S. A. and Perlman, M. D. (1998) Normal Linear Regression Models With Recursive Graphical Markov Structure. \emph{Journal of Multivariate Analysis} \bold{66}, 133--187}
\author{Gruber, Lutz F.}
\note{}


\examples{
## The predictor variable
p <- 1:10

## Two response variables
x1 <- 2*p + rnorm(10, 0, 0.1)
x2 <- cbind(50 - 2*p + rnorm(10, 0, 0.1), 0.5*p + rnorm(10, 0, 0.1))
x3 <- -p + rnorm(10, 0, 0.1)

## The regression formulas of the response variables
formulas <- list(x1 = (x1 ~ p - 1), x2 = (x2 ~ p), x3 = (x3 ~ p))

## The independence graph of the normal linear ADG model
D <- DAG(x2 ~ x1, x3 ~ x2)

## Fit the model
model <- nlADG(formulas, D)

## Check for the D-linear subspace characteristics
assertPDI(model$sigma, model$ADG, model$I) #will work
\dontrun{assertPDI(matrix(rnorm(16), ncol=4), model$ADG, model$I)} #will most likely not work
}

\seealso{\code{\link{nlADG}}, \code{\link{DAG}}, \code{\link{topSort}}}
