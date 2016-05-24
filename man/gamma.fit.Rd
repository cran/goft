\name{gamma_fit}
\alias{gamma_fit}
\title{ Fitting the Gamma distribution to data}
\description{
Fits a Gamma distribution to a random sample of positive real numbers using Villasenor and Gonzalez-Estrada (2015) parameter estimators.
}
\usage{
gamma_fit(x)
}
\arguments{
  \item{x}{  a numeric data vector containing a random sample of positive real numbers.}
}
\details{
The  Gamma distribution with shape and scale parameters is considered. 
The scale parameter is estimated by the unbiased sample estimator of the covariance  of X and log(X).
The  shape parameter is estimated by the ratio of the sample mean of X and the scale parameter estimator.
}
\value{
  Shape and scale parameter estimates.
}
\references{ Villasenor, J.A. and Gonzalez-Estrada, E. (2015). A variance ratio test of fit for Gamma distributions. \emph{Statistics and Probability Letters}, \bold{96} 1, 281-286. \url{http://dx.doi.org/10.1016/j.spl.2014.10.001}}
\author{ Elizabeth Gonzalez-Estrada \email{egonzalez@colpos.mx}, Jose A. Villasenor-Alva}
\seealso{ \code{\link{gamma_test}} for testing the Gamma distribution hypothesis.}
\examples{
# Fitting a gamma distribution to the logarithm of variable Loss contained in
# the danishuni data set 
library(fitdistrplus)
data(danishuni) 
logLoss <- log(danishuni$Loss)   # logarithm of Loss variable
logLoss <- logLoss[logLoss > 0]  # observations > 0
gamma_fit(logLoss)                 
}

