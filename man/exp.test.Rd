\name{exp_test}
\alias{exp_test}

\title{
Tests for exponentiality
}
\description{
Test based on a data trasformation and Cox-Oakes test for exponentiality.
}
\usage{
exp_test(x, method = "transf", N = 10^3)
}
\arguments{
  \item{x}{ 
 a numeric data vector containing a random sample of positive real numbers.
}
  \item{method}{ 
  test to be used. Tow available options are  \verb{"transf"} and \verb{"ratio"}. Default option is \verb{"transf"}.
}
\item{N}{ number of Monte Carlo samples used to approximate p-values when \verb{"transf"} option is chosen and the sample size is less than 200. Default is N = 10^3.
}
}
\details{
Option \verb{"transf"} performs a test based on a transformation to approximately uniformly distributed  variables.
If the sample size is larger than 200, the asymptotic null distribution of the test statistic is used to approximate the p-value; otherwise, it is approximated by Monte Carlo simulation.

Option \verb{"ratio"} performs Cox and Oakes (1984) test.
}

\value{
A list with class \code{"htest"} containing the following components.
  \item{statistic}{ the calculated value of the test statistic.  }
  \item{p.value}{   approximated p-value of the test.  }
  \item{method}{ a character string giving the name of the method used for testing exponentiality.}
  \item{data.name}{ a character string giving the name of the data set.  }
}

\references{
Cox, D.R., Oakes, D., 1984. Analysis of Survival Data. Chapman and Hall/CRC.
}
\author{
Elizabeth Gonzalez-Estrada, Jose A. Villasenor
}

\examples{
# Testing exponentiality on a simulated random sample from the exponential distribution
x <- rexp(20)
exp_test(x)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}

