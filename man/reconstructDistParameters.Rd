\name{reconstructDistParameters}
\alias{reconstructDistParameters}
\title{Reconstruct the Normal Linear ADG Model's Parameters from their D-Parameters}
\description{\code{reconstructDistParameters} reconstructs the normal linear ADG model's parameters from their \code{D}-parameters.}
\usage{
reconstructDistParameters(beta, mu, lambda, D, I, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{beta}{The \code{D}-parameter \code{beta}.}
  \item{mu}{The \code{D}-parameter \code{mu}.}
  \item{lambda}{The \code{D}-parameter \code{lambda}.}
  \item{D}{The adjacency matrix of the normal linear ADG model's independence graph sorted according to its topological order.}
  \item{I}{A list which specifies the position of the response variables in the covariance matrix.}
  \item{...}{Additional parameters.}
}
\value{
  \item{xi_hat}{The reconstructed distribution parameter \code{xi}.}
  \item{sigma_hat}{The reconstructed distribution parameter \code{sigma}.}
}
\references{Andersson, S. A. and Perlman, M. D. (1998) Normal Linear Regression Models With Recursive Graphical Markov Structure. \emph{Journal of Multivariate Analysis} \bold{66}, 133--187}
\author{Gruber, Lutz F.}
\note{The function is called by \code{nlADG}. It is not intended to be directly used by the end user.}

\seealso{\code{\link{nlADG}}, \code{\link{topSort}}}