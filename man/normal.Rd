\name{normal.test}
\alias{normal.test}
\title{ Correlation test for univariate normality }
\description{
This function implements the correlation test for normality  based on the Levy characterization by Villasenor-Alva and Gonzalez-Estrada (2013).
}
\usage{
normal.test(x)
}
\arguments{
  \item{x}{ Numeric data vector containing a random sample of  size n. }
}
\details{
Sample size (n) must be an integer between 10 and 400.

The Levy characterization of the normal distribution 
states that sums of independent normal random variables are also normal random variables.
Based on this property, the normality assumption might be judged by comparing sums of pairs of observations from a random sample to quantiles of the normal distribution. If normality holds then the pairs of sums and quantiles should lie on a straight line approximately. A formal test for normality is obtained when such comparison is based on the sample correlation coefficient of sums and quantiles.

}
\value{
  A list with class \code{"htest"} containing the following components.
  \item{statistic}{ the calculated value of the test statistic.  }
  \item{p.value}{  the p-value of the test.  }
  \item{method}{ the character string "Correlation test for normality".  }
  \item{data.name}{ a character string giving the name of the data set.  }
}
\references{ Villasenor-Alva, J.A. and Gonzalez-Estrada, E. (2013). A correlation test for normality based on the Levy characterization. \emph{Communications in Statistics: Simulation and Computation} (to appear).}

\seealso{  \code{\link{shapiro.test}}}
\examples{
x <- rnorm(20)        # generating a pseudo random sample of size 20 from
                      # the standard normal distribution
normal.test(x)        # testing normality on x
}
\keyword{ htest }



