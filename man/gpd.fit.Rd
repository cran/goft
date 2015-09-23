\name{gpd.fit}
\alias{gpd.fit}
\title{ Fitting the generalized Pareto distribution to data}
\description{
Fits a generalized Pareto distribution (gPd) to a random sample  using either the asymptotic maximum likelihood method (amle)  or the combined  estimation method (Villasenor-Alva and Gonzalez-Estrada, 2009).
}
\usage{
gpd.fit(x,method)
}
\arguments{
  \item{x}{ Numeric data vector containing a random sample of positive real number.}
  \item{method}{ a character string giving the name of the parameter estimation method to be used. There are two available methods:  \code{"combined"} and \code{"amle"}. Use \code{"combined"} for fitting a gPd with shape parameter <0. Use \code{"amle"}
                for fitting a gPd with shape parameter >= 0.  }
}
\details{
The distribution function  of the gPd is given in the details section  of the  \code{\link{gpd.test} function}.
}
\value{
  Parameter estimates.
}
\references{ Villasenor-Alva, J.A. and Gonzalez-Estrada, E. (2009). A bootstrap goodness of fit test for the generalized Pareto distribution. \emph{Computational Statistics and Data Analysis},\bold{53},11,3835-3841. \url{http://dx.doi.org/10.1016/j.csda.2009.04.001}}
\author{ Elizabeth Gonzalez-Estrada, Jose A. Villasenor-Alva }
\seealso{ \code{\link{gpd.test}} for testing the gPd hypothesis, \code{\link{rgp}} for generating gPd random numbers.}
\examples{
x <-  rgp(20,shape = 1)   	## pseudo random sample of size 20
gpd.fit(x,"amle")         ## fitting a gPd to  x  using the "amle" method
}
\keyword{ htest }
\concept{ fitting the generalized Pareto distribution to data }
