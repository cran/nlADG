\name{estimateDParameters}
\alias{estimateDParameters}
\title{Estimate the D-Parameters of a Normal Linear ADG Model}
\description{\code{estimateDParameters} calculates the MLEs of the \code{D}-parameters of a normal linear ADG model.}
\usage{
estimateDParameters(formulas, D, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formulas}{A list of regression formulas. There must be a regression formula for each response variable.}
  \item{D}{The adjacency matrix of the normal linear ADG model's independence graph sorted according to its topological order.}
  \item{...}{Additional parameters.}
}
\value{
  \item{beta_hat}{The estimate of the \code{D}-parameter \code{beta}.}
  \item{mu_hat}{The estimate of the \code{D}-parameter \code{mu}.}
  \item{lambda_hat}{The estimate of the \code{D}-parameter \code{lambda}.}
  \item{residuals}{The residuals.}
  \item{coefficients}{The estimates of the regression coefficient matrices.}
  \item{I}{A list which specifies the position of the response variables in the covariance matrix.}
  \item{models}{A list of MANOVA models corresponding to the response variables.}
  \item{fittedvalues}{The fitted values.}
}
\references{Andersson, S. A. and Perlman, M. D. (1998) Normal Linear Regression Models With Recursive Graphical Markov Structure. \emph{Journal of Multivariate Analysis} \bold{66}, 133--187}
\author{Gruber, Lutz F.}
\note{The function is called by \code{nlADG}. It is not intended to be directly used by the end user.}

\seealso{\code{\link{nlADG}}, \code{\link{topSort}}}