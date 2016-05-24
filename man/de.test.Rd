\name{laplace.test}
\alias{laplace.test}

\title{
Tests for the Laplace or double exponential distribution
}
\description{
Transformation and ratio tests for the Laplace distribution by  Gonzalez-Estrada and Villasenor (2016).
}
\usage{
laplace.test(x, method = "transf", N = 10^5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a numeric data vector containing a random sample of real numbers.
}
  \item{method}{ the type of test to be performed. Two available options are \verb{"transf"} and
  \verb{"ratio"}. Default option is \verb{"transf"}. 
  }
  \item{N}{
number of Monte Carlo samples used to approximate the p-value of the test when the \verb{"ratio"} option is chosen and the sample size is less than 500. Default is N = 10^5.
}
}
\details{

When \verb{"transf"} option is chosen, a transformation to approximately exponential random variables is performed and the exponentiality hypothesis is assessed using Anderson-Darling  test. 

When \verb{"ratio"} option is chosen, a test based on the ratio of two estimators of the scale parameter is performed. For samples of size n < 500, the p-value of this test is approximated by Monte Carlo simulation. Otherwise, it is approximated by the standard normal cumulative distribution function. 
}
\value{
  A list with class \code{"htest"} containing the following components.
  \item{statistic}{ the calculated value of the test statistic.  }
  \item{p.value}{   approximated p-value of the test.  }
  \item{method}{ a character string giving the name of the method used for testing the null hypothesis.}
  \item{data.name}{ a character string giving the name of the data set.  }
}
\references{
Gonzalez-Estrada, E. and Villasenor, J.A. (2016). A ratio goodness-of-fit test for the Laplace distribution. (Submitted).
}
\author{
Elizabeth Gonzalez-Estrada \email{egonzalez@colpos.mx}, Jose A. Villasenor
}

\examples{
# Example 1:  testing the Laplace distribution hypothesis using "transf" option
x <- rnorm(50)   # simulating a random sample from a normal distribution
laplace.test(x)    

# Example 2: testing the Laplace distribution hypothesis using "ratio" option
x <- rt(60,4)    # simulating a random sample from Student's t distribution with 4 d.f.
laplace.test(x, method = "ratio")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }

