\name{mvShapiro.Test}
\alias{mvShapiro.Test}
\title{ Shapiro-Wilk test for multivariate normality}
\description{
Villasenor-Alva and Gonzalez-Estrada (2009) generalization of Shapiro-Wilk test for multivariate normality.
}
\usage{
mvShapiro.Test(X)
}
\arguments{
  \item{X}{ numeric data matrix with d columns (vector dimension) and n rows (sample size). }
}
\details{
Sample size (n) must be larger than the number of coordinates (d).

When d=1, \code{mvShapiro.Test(X)} produces the same results as \code{shapiro.test(X)}.
}
\value{
  A list with class \code{"htest"} containing the following components.
  \item{statistic}{ the value of the generalized Shapiro-Wilk statistic for testing multivariate normality.  }
  \item{p.value}{  the approximated  p-value of the test.  }
  \item{method}{ the character string "Generalized Shapiro-Wilk test for multivariate normality".  }
  \item{data.name}{ a character string giving the name of the data set.  }
}
\references{ Villasenor-Alva, J.A. and Gonzalez-Estrada, E. (2009). A generalization of Shapiro-Wilk's test for multivariate normality. \emph{Communications in Statistics: Theory and Methods}, \bold{38} 11, 1870-1883. \url{http://dx.doi.org/10.1080/03610920802474465}}
\author{ Elizabeth Gonzalez-Estrada \email{egonzalez@colpos.mx}, Jose A. Villasenor-Alva}
\seealso{ \code{\link{shapiro.test}} and \code{\link{normal.test}}  for testing univariate normality.}
\examples{
#-----------------------------------------------------------------------------------
# iris.virginica contains a set of measurements corresponding to 
# Iris virginica of the famous  iris data set.

iris.virginica <- as.matrix(iris[iris$Species == "virginica",1:4],ncol=4) 
mvShapiro.Test(iris.virginica)    # Testing multivariate normality on iris.virginica

#-----------------------------------------------------------------------------------
# Testing multivariate normality on the goats dataset.
data(goats)
mvShapiro.Test(as.matrix(goats))
}
\keyword{ htest }
\keyword{ multivariate analysis}