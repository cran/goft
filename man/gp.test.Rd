\name{gp.test}
\alias{gp.test}
\title{ Bootstrap test for the generalized Pareto distribution}
\description{
Test of fit for the generalized Pareto distribution (gPd) with unknown parameters by Villasenor-Alva and Gonzalez-Estrada (2009).
}
\usage{
gp.test(x, B = 999)
}
\arguments{
  \item{x}{ numeric data vector containing a random sample of positive real numbers.}
  \item{B}{ number of bootstrap samples used to approximate p-values. Default is \code{B=999}.}
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
A list with class \code{"htest"} containing the following components.
  \item{p.value}{  an approximated p-value of the test using parametric bootstrap.  }
  \item{method}{ the character string "Bootstrap test of fit for the generalized Pareto distribution".  }
  \item{data.name}{ a character string giving the name of the data set.  }
  \item{pvalues}{ approximated p-values of the tests for \eqn{H_0^-}{H_0^-} and \eqn{H_0^+}{H_0^+}}
}
\references{ Villasenor-Alva, J.A. and Gonzalez-Estrada, E. (2009). A bootstrap goodness of fit test for the generalized Pareto distribution. \emph{Computational Statistics and Data Analysis},\bold{53},11,3835-3841.  \url{http://dx.doi.org/10.1016/j.csda.2009.04.001}}
\author{ Elizabeth Gonzalez-Estrada \email{egonzalez@colpos.mx}, Jose A. Villasenor-Alva}
\seealso{  \code{\link{gp.fit}} for fitting a gPd to data.}
\examples{
# Testing the gPd hypothesis on the excesses above the threshold 0.165 ppm of the ozone
# levels given in the o3 data set
data(o3)
o3levels <- o3$ozone_level - 0.165  # ozone levels minus the threshold 0.165 ppm     
gp.test(o3levels, B = 9999)   # testing the gPd hypothesis using 9999 bootstrap samples            
}

\keyword{htest}

