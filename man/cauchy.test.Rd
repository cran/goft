\name{cauchy_test}
\alias{cauchy_test}

\title{
Tests for the Cauchy distribution
}
\description{
Two tests for the Cauchy distribution hypothesis.
}
\usage{
cauchy_test(x, N = 10^3, method = "transf")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a numeric data vector containing a random sample of real numbers.
}
  \item{N}{
number of Monte Carlo samples used to approximate the p-value of the tests. Default is N = 10^3.
}
\item{method}{ character string giving the name of the method to be used for testing the Cauchy distribution hypothesis. Two available options are \verb{"transf"} and \verb{"ratio"}.
}
}
\details{
Option  \verb{"ratio"} performs a test for the Cauchy distribution based on the ratio of the maximum likelihood estimator for the scale parameter  and the mean absolute deviation (Gonzalez-Estrada and Villasenor, 2018). 
 
  Option \verb{"transf"} performs a test based on a data transformation to    approximately exponentially distributed data  (Villasenor and Gonzalez-Estrada, 2020).
 
}
\value{
  A list with class \code{"htest"} containing the following components.
  \item{statistic}{ the calculated value of the test statistic.  }
  \item{p.value}{  the approximated p-value of the test.  }
  \item{method}{ the character string "Test for the Cauchy distribution based on the ratio of two scale estimators".}
  \item{data.name}{ a character string giving the name of the data set.  }
}

\references{ Gonzalez-Estrada, E. and Villasenor, J.A. (2018). An R package for testing goodness of fit: goft. \emph{Journal of Statistical Computation and Simulation}, \bold{88} 4, 726-751. \url{https://doi.org/10.1080/00949655.2017.1404604}

Villasenor, J.A. and Gonzalez-Estrada, E.  (2020). Goodness of fit tests for Cauchy distributions using data transformations. In I. Ghosh, N. Balakrishnan and H.K.T. Ng. \emph{Contributions of Barry C. Arnold to Statistical Science - Theory and Applications}. Springer. 

}
\author{ Elizabeth Gonzalez-Estrada \email{egonzalez@colpos.mx}}

\examples{
x <- rnorm(20)    # simulating a data set from a normal distribution
cauchy_test(x)    # testing the Cauchy distribution hypothesis 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }

