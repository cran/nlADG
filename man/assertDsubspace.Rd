\name{assertDsubspace}
\alias{assertDsubspace}
\title{Assert that the D-Linear Subspace Properties Hold}
\description{\code{assertDsubspace} compares the design matrices of the given formulas to hint for a possible violation of the \code{D}-linear subspace characteristics.}
\usage{
assertDsubspace(formulas, D, ...)
}
\arguments{
  \item{formulas}{A list of regression formulas. There must be a regression formula for each response variable.}
  \item{D}{The adjacency matrix of the normal linear ADG model's independence graph.}
  \item{...}{Additional parameters.}
}
\details{The function \code{assertDsubspace} produces an error if the regression system does not admit the \code{D}-linear subspace characteristics.}
\references{Andersson, S. A. and Perlman, M. D. (1998) Normal Linear Regression Models With Recursive Graphical Markov Structure. \emph{Journal of Multivariate Analysis} \bold{66}, 133--187}
\author{Gruber, Lutz F.}
\note{}


\examples{
## The predictor variable
p <- 1:10

## Two response variables
x1 <- 2*p + rnorm(10, 0, 0.1)
x2 <- cbind(50 - 2*p + rnorm(10, 0, 0.1), 0.5*p + rnorm(10, 0, 0.1))

## The regression formulas of the response variables
formulas <- list(x1 = (x1 ~ p - 1), x2 = (x2 ~ p))

## Check for the D-linear subspace characteristics
assertDsubspace(formulas, D = DAG(x2 ~ x1)) #will go through
\dontrun{assertDsubspace(formulas, D = DAG(x1 ~ x2))} #will produce an error
}

\seealso{\code{\link{nlADG}}, \code{\link{DAG}}}
