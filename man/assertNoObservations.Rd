\name{assertNoObservations}
\alias{assertNoObservations}
\title{Assert that Enough Observations Are Given}
\description{\code{assertNoObservations} checks if enough observations are given for the maximum likelihood estimates to exist and be unique. If not, the function produces an error.}
\usage{
assertNoObservations(formulas, D, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formulas}{A list of regression formulas. There must be a regression formula for each response variable.}
  \item{D}{The adjacency matrix of the normal linear ADG model's independence graph sorted according to its topological order.}
  \item{...}{Additional parameters.}
}
\references{Andersson, S. A. and Perlman, M. D. (1998) Normal Linear Regression Models With Recursive Graphical Markov Structure. \emph{Journal of Multivariate Analysis} \bold{66}, 133--187}
\author{Gruber, Lutz F.}
\note{}


\examples{
## The predictor variable
v <- 1:5

## The response variable
x <- 25 - 2*v + rnorm(5)

## The regression formula of the response variable
formulas <- list(x=(x~v))

## The model's independence graph
D <- DAG(x~x)

## Check if enough observations are given
assertNoObservations(formulas, D) #will work

## Reduce the number of observations
x <- x[1:2]
v <- v[1:2]
\dontrun{assertNoObservations(formulas, D)} #will produce an error message
}

\seealso{\code{\link{nlADG}}, \code{\link{DAG}}, \code{\link{topSort}}}
