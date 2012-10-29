\name{selectInd}
\alias{selectInd}
\title{Select Entries of the MLE sigma Matrix by their Names}
\description{Given the index list \code{I}, the \code{selectInd} function allows for the convenient selection of entries in the \code{sigma} matrix by the names of the respective response variables.}
\usage{
selectInd(indexList, selectedIndices)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{indexList}{A list that specifies the rows and columns in which the values to the response variables can be found.}
  \item{selectedIndices}{A vector with the names of the response variables whose entries to select.}
}
\value{The row and column indices of the selected response variables in the \code{sigma} matrix.}
\references{Andersson, S. A. and Perlman, M. D. (1998) Normal Linear Regression Models With Recursive Graphical Markov Structure. \emph{Journal of Multivariate Analysis} \bold{66}, 133--187}
\author{Gruber, Lutz F.}
\note{}


\examples{
## Set up a normal linear ADG model
v <- 1:100
x <- 2*v + rnorm(100, 0, 5)
y <- x + 50 + rnorm(100, 0, 1)

D <- DAG(y~x)
formulas <- list(x=(x~v-1), y=(y~v))

model <- nlADG(formulas, D)

## The MLE of sigma
sigma <- model$sigma

## The indices of the response variable x
xInd <- selectInd(model$I, "x")
sigma[xInd, xInd]

## The indices to the response variables pa("y", model$ADG)
ypInd <- selectInd(model$I, pa("y", model$ADG))
sigma[ypInd, ypInd]

## The indices to the response variables c("x", "y")
xypInd <- selectInd(model$I, c("x", "y"))
sigma[xypInd, xypInd]
}

\seealso{\code{\link{nlADG}}, \code{\link{DAG}}, \code{\link{pa}}}