\name{weibull_test}
\alias{weibull_test}

\title{
Test for the Weibull distribution 
}
\description{
Test for the Weibull distribution based on a data transformation to observations with  type I extreme value (Gumbel) distribution.
}
\usage{
weibull_test(x, method = "transf", N = 1000)
}

\arguments{
\item{x}{ numeric data vector containing a random sample of positive observations.}
  \item{method }{ test to be used. So far,   \verb{"transf"}  is the only option and it comes by default.}
  \item{N}{ number of Monte Carlo samples used to approximate the p-value of the test. Default is N = 1000.}
}
\details{
The \verb{"ratio"} test provided by function \code{\link{ev_test}} for the type I extreme value distribution  is applied to the transformed observations.
}
\value{
A list with class \code{"htest"} containing the following components.
  \item{statistic}{ the calculated value of the test statistic.  }
  \item{p.value}{  an approximated p-value of the test.  }
  \item{method}{ the character string "Test for the Weibull distribution".  }
  \item{data.name}{ a character string giving the name of the data set.  }
}

\author{
Elizabeth Gonzalez-Estrada, Jose A. Villasenor 
}

\seealso{
 Other tests for positive datasets are provided by functions \code{\link{ig_test}}, \code{\link{gamma_test}} and \code{\link{lnorm_test}}.
}
\examples{
# Testing the Weibull distribution hypothesis on the maximum ozone levels given
# in the o3max data set.
data(o3max)
x <- o3max[ ,2]     # maximum ozone levels
weibull_test(x, N = 10000)     
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }
 