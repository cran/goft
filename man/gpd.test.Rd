\name{gpd.test}
\alias{gpd.test}
\title{ Bootstrap test for the generalized Pareto distribution}
\description{
Villasenor-Alva and Gonzalez-Estrada (2009) test of fit for the generalized Pareto distribution (gPd) with unknown parameters.
}
\usage{
gpd.test(x,J)
}
\arguments{
  \item{x}{ numeric data vector containing a random sample of positive real numbers.}
  \item{J}{ number of bootstrap samples. This is an optional argument. Default \code{J=999}.}
}
\details{
This bootstrap test for the null hypothesis \eqn{H_0:}{H_0:} a random sample has a  gPd with unknown shape parameter \eqn{\gamma}{gamma} is an intersection-union test for the hypotheses \eqn{H_0^-:}{H_0^-:} a random sample has a gPd with  \eqn{\gamma < 0}{gamma <0 }, and \eqn{H_0^+:}{H_0^+:} a random sample has a gPd with  \eqn{\gamma >=0}{gamma >=0}.
Thus, heavy and non-heavy tailed gPd's are included in the null hypothesis. The parametric bootstrap is performed on \eqn{\gamma}{gamma} for each of the two hypotheses.

The gPd function with unknown shape and scale parameters  \eqn{\gamma}{gamma} and \eqn{\sigma}{sigma}  is given by

  \deqn{F(x) = 1 - \left[ 1 + \frac{\gamma x}{ \sigma } \right] ^ { - 1 /\gamma},}{ F(x) = 1 - [ 1 + gamma x  / sigma ]^(-1/gamma),}

  where   \eqn{\gamma}{gamma}  is a real number, \eqn{\sigma > 0}{sigma > 0} and \eqn{1 + \gamma  x  / \sigma > 0}{1 + gamma x / sigma > 0}. When \eqn{\gamma = 0}{gamma =
    0}, F(x) becomes the exponential distribution with scale parameter \eqn{\sigma}{sigma}:   \deqn{F(x) = 1 -exp\left(-x/\sigma \right).}{1-exp(-x/sigma).}
}
\value{
  A list with the following components.
  \item{boot.test}{  a list with class \code{"htest"}  containing the  p-value of the test, the name of the data set,  and  the character string "Bootstrap test of fit for the generalized Pareto distribution".  }
  \item{p.values}{ the p-values of the tests of the hypotheses \eqn{H_0^-}{H_0^-} and \eqn{H_0^+}{H_0^+} described above. }
}
\references{ Villasenor-Alva, J.A. and Gonzalez-Estrada, E. (2009). A bootstrap goodness of fit test for the generalized Pareto distribution. \emph{Computational Statistics and Data Analysis},\bold{53},11,3835-3841.  \url{http://dx.doi.org/10.1016/j.csda.2009.04.001}}
\author{ Elizabeth Gonzalez-Estrada \email{egonzalez@colpos.mx}, Jose A. Villasenor-Alva}
\seealso{  \code{\link{gpd.fit}} for fitting a gPd to data,  \code{\link{rgp}} for generating  gPd random numbers.}
\examples{
x <-  rgp(20,shape = 1)   	## generating a random sample of size 20 from a gPd
gpd.test(x)                 ## testing the gPd hypothesis on x
}

\concept{ hypothesis testing}
\concept{ extreme value analysis }
\concept{ parametric bootstrap}
