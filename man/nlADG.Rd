\name{nlADG}
\alias{nlADG}
\alias{nlADG.default}
\alias{print.nlADG}
\alias{summary.nlADG}
\alias{print.summary.nlADG}
\alias{predict.nlADG}
\title{Regression in the Normal Linear ADG Model}
\description{\code{nlADG} calculates the MLEs of the normal linear ADG model's parameters.}
\usage{
nlADG(formulas, ...)

\method{nlADG}{default}(formulas, D, ...)
\method{print}{nlADG}(x, ...)
\method{summary}{nlADG}(object, ...)
\method{print.summary}{nlADG}(x, ...)
\method{predict}{nlADG}(object, newdata, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formulas}{A list of regression formulas. There must be a regression formula for each response variable.}
  \item{D}{The adjacency matrix of the normal linear ADG model's independence graph.}
  \item{x}{An object of the \code{nlADG} class.}
  \item{object}{An object of the \code{nlADG} class.}
  \item{newdata}{The data for which to predict the outcome.}
  \item{...}{Additional parameters, which may be used by functions deeper in the structure of this program.}
}
\details{For any entry in \code{formulas}, the name of the entry has to match the name of the corresponding response variable.}
\value{
  \item{coefficients}{The estimates of the regression coefficients.}
  \item{sigma}{The estimate of the covariance matrix.}
  \item{xi}{The estimate of the mean of the distribution.}
  \item{fitted.values}{The fitted values given the predictors.}
  \item{residuals}{The residuals.}
  \item{MANOVAmodels}{A list of MANOVA models corresponding to the response variables.}
  \item{ADG}{The sorted adjacency matrix of the model's independence graph \code{D}.}
  \item{I}{A list which specifies the position of the response variables in the covariance matrix.}
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

## The regression formulas of the response variables
formulas <- list(x1 = (x1 ~ p - 1), x2 = (x2 ~ p))

## The normal linear ADG model's independence graph
D <- DAG(x2 ~ x1)

## Fit the model
model <- nlADG(formulas, D)


model$xi # Access the MLE of xi
model$sigma # Access the MLE of sigma

## Predict values
newdata <- data.frame(p=100:105)
predict(model, newdata)
}

\seealso{\code{\link{DAG}}}

%\keyword{ADG}
%\keyword{DAG}
%\keyword{Regression}
%\keyword{Normal Linear ADG Model}% __ONLY ONE__ keyword per line
