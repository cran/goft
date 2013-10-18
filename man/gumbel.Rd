\name{gumbel.test}
\alias{gumbel.test}

\title{ Correlation test for the Gumbel distribution }
\description{
This function implements the correlation test of fit for the Gumbel distribution (also known as  type I extreme value distribution) based on the max-stability property (Gonzalez-Estrada and Villasenor-Alva,  2010).
}
\usage{
gumbel.test(x)
}
\arguments{
  \item{x}{ Numeric data vector containing a random sample of  size n. }
}
\details{
Sample size (n)  must be an integer between 20 and 250.

The max-stability property of the Gumbel distribution states that the maximum of a number of Gumbel random variables is also a Gumbel random variable. 

The test statistic is the sample correlation of the maximum of pairs of observations and quantiles of the Gumbel distribution. 
}
\value{
  A list with class \code{"htest"} containing the following components.
  \item{statistic}{ the calculated value of the test statistic.  }
  \item{p.value}{  the p-value of the test.  }
  \item{method}{ the character string "Test of fit for the Gumbel distribution".  }
  \item{data.name}{ a character string giving the name of the data set.  }
}
\references{ Gonzalez-Estrada, E. and Villasenor-Alva, J.A. (2010). A Goodness-of-Fit Test for Location-Scale Max-Stable Distributions. \emph{Communications in Statistics: Simulation and Computation}, \bold{39} 3, 557-562.  }


\examples{
x <- -log(rexp(20))   # generating a pseudo random sample of size 20 from 
                      # the standard Gumbel distribution
gumbel.test(x)        # testing the Gumbel hypothesis on x
}
\keyword{ htest }

