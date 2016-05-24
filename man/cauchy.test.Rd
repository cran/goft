\name{cauchy.test}
\alias{cauchy.test}

\title{
Test for the Cauchy distribution
}
\description{
A test based on the ratio of two estimators for the scale parameter of the Cauchy distribution.
}
\usage{
cauchy.test(x, N = 10^3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a numeric data vector containing a random sample of real numbers.
}
  \item{N}{
number of Monte Carlo samples used to approximate the p-value of the test. Default is N = 10^3.
}
}
\details{

The test statistic is the ratio of the maximum likelihood estimator of the scale parameter and the mean absolute deviation. 

}
\value{
  A list with class \code{"htest"} containing the following components.
  \item{statistic}{ the calculated value of the test statistic.  }
  \item{p.value}{  the approximated p-value of the test.  }
  \item{method}{ the character string "Test for the Cauchy distribution based on the ratio of two scale estimators".}
  \item{data.name}{ a character string giving the name of the data set.  }
}

\author{
Elizabeth Gonzalez-Estrada and Jose A. Villasenor
}

\examples{
x <- rnorm(20)    # simulating a data set from a normal distribution
cauchy.test(x)    # testing the Cauchy distribution hypothesis 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }

