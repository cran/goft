\name{gp_fit}
\alias{gp_fit}
\title{ Fitting the generalized Pareto distribution to data}
\description{
Fits a generalized Pareto distribution (gPd) to a random sample  using either the asymptotic maximum likelihood method (amle)  or the combined  estimation method (Villasenor-Alva and Gonzalez-Estrada, 2009).
}
\usage{
gp_fit(x, method)
}
\arguments{
  \item{x}{ a numeric data vector containing a random sample of positive real numbers.}
  \item{method}{ a character string giving the name of the parameter estimation method to be used. There are two available methods:  \code{"combined"} and \code{"amle"}. Use \code{"combined"} for fitting a gPd with negative shape parameter. Use \code{"amle"} for fitting a gPd with non-negative shape parameter.  }
}
\details{
The distribution function  of the gPd is provided in the details section  of the  \code{\link{gp_test} function}.


}
\value{
  Parameter estimates.
}
\references{ Villasenor-Alva, J.A. and Gonzalez-Estrada, E. (2009). A bootstrap goodness of fit test for the generalized Pareto distribution. \emph{Computational Statistics and Data Analysis},\bold{53},11,3835-3841. \url{http://dx.doi.org/10.1016/j.csda.2009.04.001}}
\author{ Elizabeth Gonzalez-Estrada, Jose A. Villasenor-Alva }
\seealso{ \code{\link{gp_test}} for testing the gPd hypothesis.}
\examples{
# Example 1: fitting a gPd to the excesses above the threshold 0.165 ppm of the ozone
# levels provided in the o3 data set.
data(o3)
o3levels <- o3$ozone_level - 0.165      # ozone levels minus the threshold 0.165 ppm 
gp_fit(o3levels, method = "combined")   # fitting the gPd using the "combined" method 

# Example 2: fitting a gPd to a simulated data set from the uniform distribution,
# which is a gPd with shape parameter = -1 and scale parameter = 1.
x <-  runif(100)       # simulating a uniform random sample 
gp_fit(x, "combined")  # fitting a gPd to  x  using the "combined" method
}
\keyword{ htest }

