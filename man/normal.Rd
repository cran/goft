\name{normal_test}
\alias{normal_test}
\title{ Test for normality }
\description{
Correlation test of fit for normality  based on the Levy characterization (Villasenor-Alva and Gonzalez-Estrada, 2015).
}
\usage{
normal_test(x, method = "cor")
}
\arguments{
  \item{x}{ a numeric data vector containing a random sample of  size n. }
  \item{method}{ a character string giving the name of the test to be used. So far the only option is \verb{"cor"}. }
}
\details{
Sample size (n) must be an integer lying between 10 and 400.

The Levy characterization of the normal distribution 
states that sums of independent normal random variables are also normal random variables.
Based on this property, the normality assumption might be judged by comparing sums of pairs of observations from a random sample to quantiles of the standard normal distribution. If normality holds then the pairs of sums and quantiles should lie on a straight line approximately. A formal test for normality is obtained when such a comparison is based on the sample correlation coefficient of sums and quantiles.

}
\value{
  A list with class \code{"htest"} containing the following components.
  \item{statistic}{ the calculated value of the test statistic.  }
  \item{p.value}{  an approximated p-value of the test.  }
  \item{method}{ the character string "Correlation test for normality".  }
  \item{data.name}{ a character string giving the name of the data set.  }
}
\references{ Villasenor-Alva, J.A. and Gonzalez-Estrada, E. (2015). A correlation test for normality based on the Levy characterization. \emph{Communications in Statistics: Simulation and Computation}, \bold{44} 5, 1225-1238. \url{http://dx.doi.org/10.1080/03610918.2013.810261}}
\author{ Elizabeth Gonzalez-Estrada \email{egonzalez@colpos.mx}, Jose A. Villasenor-Alva}
\seealso{  \code{\link{shapiro.test}}}
\examples{
data(goats)                   # loading the "goats" data set
apply(goats,2,normal_test)    # testing normality on each variable of the "goats" data set    
}
\keyword{ htest }



