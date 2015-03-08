\name{gammadist.test}
\alias{gammadist.test}
\title{ Test for the gamma distribution }
\description{
Test of fit for the gamma distribution proposed by Villasenor and Gonzalez-Estrada (2015) \url{http://dx.doi.org/10.1016/j.spl.2014.10.001}.
}
\usage{
gammadist.test(x)
}
\arguments{
  \item{x}{ Numeric data vector containing a random sample of  size n. }
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
  \item{method}{ the character string "Test for the gamma distribution".  }
  \item{data.name}{ a character string giving the name of the data set.  }
}
\references{ Villasenor, J.A. and Gonzalez-Estrada, E. (2015). A variance ratio test of fit for gamma distributions. \emph{Statistics and Probability Letters}, \bold{96} 1, 281-286.\url{http://dx.doi.org/10.1016/j.spl.2014.10.001}}
\author{ Elizabeth Gonzalez-Estrada \email{egonzalez@colpos.mx}, Jose A. Villasenor-Alva }
\examples{
x <- rgamma(50,3)
gammadist.test(x)        # testing the gamma hypothesis
}




