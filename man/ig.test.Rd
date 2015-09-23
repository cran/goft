\name{ig.test}
\alias{ig.test}
\title{ A transformation test for Inverse Gaussian distributions}
\description{

Test of fit for  Inverse Gaussian distributions with unknown location and scale parameters based on a transformation to gamma variables (Villasenor and Gonzalez-Estrada, 2015).

}
\usage{
ig.test(x)
}
\arguments{
  \item{x}{ numeric data vector containing a random sample of  positive real numbers. }
}
\details{
Anderson-Darling test for a transformation of data to gamma variables with shape parameter equal to 1/2  is used.
 
The right-hand side tail of the null distribution for Anderson-Darling test statistic is approximated by an inverse Gaussian distribution in order to obtain approximated p-values. 

}
\value{
  A list with class \code{"htest"} containing the following components.
  \item{statistic}{ the calculated value of the test statistic.  }
  \item{p.value}{  the approximated p-value of the test.  }
  \item{method}{ the character string "A transformation test for Inverse Gaussian distributions".}
  \item{data.name}{ a character string giving the name of the data set.  }
}
\references{ Villasenor, J.A. and Gonzalez-Estrada, E. (2015). Tests of fit for Inverse Gaussian distributions. \emph{Statistics and Probability Letters}, \bold{105}, 189-194. \url{http://dx.doi.org/10.1016/j.spl.2015.06.017}}
\author{ Elizabeth Gonzalez-Estrada \email{egonzalez@colpos.mx}, Jose A. Villasenor-Alva}
\examples{
x <- rgamma(50,2)
ig.test(x)        # testing the composite Inverse Gaussian hypothesis 
}

\keyword{ htest }


