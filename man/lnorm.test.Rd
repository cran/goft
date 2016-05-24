\name{lnorm.test}
\alias{lnorm.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Test for the lognormal distribution
}
\description{
Test for the lognormal distribution based on a data transformation to normal observations.
}
\usage{
lnorm.test(x)
}
\arguments{
  \item{x}{ a numeric data vector containing a random sample of positive observations.}
}
\details{
Shapiro-Wilk test is used for testing normality on the transformed observations.
}
\value{
  A list with class \code{"htest"} containing the following components.
  \item{statistic}{ the calculated value of the test statistic.  }
  \item{p.value}{  an approximated p-value of the test.  }
  \item{method}{ the character string "Test for the lognormal distribution based on a transformation to normality".  }
  \item{data.name}{ a character string giving the name of the data set.  }
}

\author{
Elizabeth Gonzalez-Estrada, Jose A. Villasenor 
}

\seealso{
 Other tests for some alternative distributions to the lognormal model are impremented in  functions \code{\link{ig.test}}, \code{\link{gamma_test}} and \code{\link{weibull.test}}.
}
\examples{
# Testing the lognormal distribution hypothesis  on the compressive strength variable
# of the strength data set.
data("strength")
x  <- strength$cstrength   # compressive strength
lnorm.test(x)              # testing the lognormal distribution hypothesis 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }

