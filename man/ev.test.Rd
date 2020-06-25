\name{ev_test}
\alias{ev_test}

\title{
Tests for the extreme value distributions
}
\description{
Correlation and ratio tests for Fisher-Tippet extreme value distributions.
}
\usage{
ev_test(x, dist = "gumbel", method = "cor", N = 1000)
}

\arguments{
  \item{x}{ a numeric data vector containing a random sample. }
  \item{dist}{ the extreme value distribution to be tested. Options are \verb{"gumbel"}, \verb{"frechet"} and \verb{"weibull"} for testing the Gumbel, Frechet and Weibull distribution hypotheses. Default option is \verb{"gumbel"}.} 
  \item{method}{ the test to be used. Available options are \verb{"ratio"} and \verb{"cor"}. Default is \verb{"cor"}.} 
  \item{N}{ number of Monte Carlo samples used to approximate the p-value of the test when \verb{"ratio"} option is chosen. Default is N = 1000. }
}

\details{
Option  \verb{"ratio"} performs a test based on the ratio of two estimators for the variance of the Gumbel (type I extreme value) distribution (Gonzalez-Estrada and Villasenor, 2018).

Option \verb{"cor"} performs a test based on the max-stability property of extreme value distributions (Gonzalez-Estrada and Villasenor, 2010). The sample size must lie between 20 and 250.
}
\value{
  A list with class \code{"htest"} containing the following components.
  \item{statistic}{ the calculated value of the test statistic.  }
  \item{p.value}{  the approximated p-value of the test.  }
  \item{method}{ a character string for the method used to test the null hypothesis.  }
  \item{data.name}{ a character string giving the name of the data set.  }
}
\references{
Gonzalez-Estrada, E. and Villasenor, J.A. (2018). An R package for testing goodness of fit: goft. \emph{Journal of Statistical Computation and Simulation}, \bold{88} 4, 726-751. \url{https://doi.org/10.1080/00949655.2017.1404604}

 Gonzalez-Estrada, E. and Villasenor-Alva, J.A. (2010). A Goodness-of-Fit Test for Location-Scale Max-Stable Distributions. \emph{Communications in Statistics: Simulation and Computation}, \bold{39} 3, 557-562. 
 \url{https://doi.org/10.1080/03610910903528293}
 
}
\author{
Elizabeth Gonzalez-Estrada, Jose A. Villasenor
}


\examples{
# Example 1:  testing the Gumbel distribution hypothesis using the correlation test
x <- -log(rexp(20))   # simulating a data set from the standard Gumbel distribution
ev_test(x)            

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }

