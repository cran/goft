\name{gammadist.test}
\alias{gammadist.test}
\title{ Test for the Gamma distribution }
\description{
Villasenor and Gonzalez-Estrada (2015) test of fit for the Gamma distribution with unknown shape and scale parameters.
}
\usage{
gammadist.test(x)
}
\arguments{
  \item{x}{ Numeric data vector containing a random sample  of positive real numbers. }
}
\details{
The test statistic is based on the ratio of two estimators of the variance, namely, the sample variance and the moments estimator obtained by Villasenor and Gonzalez-Estrada (2015), which is the product of the sample mean of X by the sample covariance of X and log(X).

The asymptotic null distribution of the test statistic is used to approximate p-values.

NOTE: the unbiased sample covariance estimator is used to compute the test statistic.

}
\value{
  A list with class \code{"htest"} containing the following components.
  \item{statistic}{ the calculated value of the test statistic.  }
  \item{p.value}{  the approximated p-value of the test.  }
  \item{method}{ the character string "Test of fit for the Gamma distribution".  }
  \item{data.name}{ a character string giving the name of the data set.  }
}
\references{ Villasenor, J.A. and Gonzalez-Estrada, E. (2015). A variance ratio test of fit for Gamma distributions. \emph{Statistics and Probability Letters}, \bold{96} 1, 281-286. \url{http://dx.doi.org/10.1016/j.spl.2014.10.001}}
\author{ Elizabeth Gonzalez-Estrada \email{egonzalez@colpos.mx}, Jose A. Villasenor-Alva }
\seealso{  \code{\link{gammadist.fit}} for fitting a Gamma distribution to data.}
\examples{
x <- rgamma(50,3)
gammadist.test(x)        # testing the gamma hypothesis
}




