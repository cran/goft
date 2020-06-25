\name{ig_test}
\alias{ig_test}
\title{ Tests for the Inverse Gaussian distribution}
\description{
Implementation of three tests of fit for  Inverse Gaussian distributions with unknown parameters (see Details section).

}
\usage{
ig_test(x, method = "transf")
}
\arguments{
  \item{x}{ a numeric data vector containing a random sample of  positive real numbers. }
  \item{method}{ character string giving the name of the method to be used for testing the Inverse Gaussian hypothesis. Two available options are \verb{"transf"} and \verb{"ratio"}.}
}
\details{
Option \verb{"transf"} performs two tests based on data transformations. One test 
transforms the observations to approximately normally distributed observations and then uses Shapiro-Wilk test for asssessing univariate normality  (Villasenor et al., 2019). 
The other test is based on a transformation of data to gamma variables with shape parameter equal to 1/2 and uses Anderson-Darling test for testing the gamma distribution (Villasenor and Gonzalez-Estrada, 2015).
The right-hand side tail of the null distribution for Anderson-Darling test statistic is approximated by an inverse Gaussian distribution in order to get approximated p-values. 

Option \verb{"ratio"} performs a test based on the ratio of two variance estimators. It is based on the asymptotic distribution of the test statistic. The use of this test is appropriate for large sample sizes (Villasenor and Gonzalez-Estrada, 2015).
}
\value{
  A list with class \code{"htest"} containing the following components.
  \item{statistic}{ the calculated value of the test statistic.  }
  \item{p.value}{  the approximated p-value of the test.  }
  \item{method}{ a character string giving the name of the method used for testing the null hypothesis.}
  \item{data.name}{ a character string giving the name of the data set.  }
}
\references{ Villasenor, J.A. and Gonzalez-Estrada, E. (2015). Tests of fit for Inverse Gaussian distributions. \emph{Statistics and Probability Letters}, \bold{105}, 189-194. \url{http://dx.doi.org/10.1016/j.spl.2015.06.017}

 Villasenor, J.A. Gonzalez-Estrada, E. and  Ochoa, A. (2019). On Testing the Inverse Gaussian Distribution Hypothesis. \emph{ Sankhya B:  The Indian Journal of Statistics}. \bold{81}, 60-74. \url{https://doi.org/10.1007/s13571-017-0148-8}
}

\author{ Elizabeth Gonzalez-Estrada \email{egonzalez@colpos.mx}, Jose A. Villasenor-Alva}
\examples{
# Example 1: testing the inverse Gaussian distribution hypothesis  using "transf" 
# option on the compressive strength variable of the strength data set.
data("strength")
y <- strength$cstrength        # compressive strength
ig_test(y, method = "transf")

# Getting the p-value of the test based on a transformation to gamma variables
ig_test(y)[[2]]$p.value 

# Example 2: testing the inverse Gaussian distribution hypothesis  using "ratio" 
# option on a simulated data set.
x <- rlnorm(500)    # simulating a random sample from the lognormal distribution
ig_test(x, method = "ratio") 

}

\keyword{ htest }


