\name{mvshapiro_test}
\alias{mvshapiro_test}
\title{ Shapiro-Wilk test for multivariate normality}
\description{
 A generalization of Shapiro-Wilk test for multivariate normality (Villasenor-Alva and Gonzalez-Estrada, 2009).
}
\usage{
mvshapiro_test(X)
}
\arguments{
  \item{X}{ a numeric data matrix with d columns  and n rows. }
}
\details{
Sample size (n) must be larger than  vector dimension (d).

When d = 1, \code{mvshapiro_test(X)} produces the same results as \code{shapiro.test(X)}.
}
\value{
  A list with class \code{"htest"} containing the following components.
  \item{statistic}{ the value of the generalized Shapiro-Wilk statistic for testing multivariate normality.  }
  \item{p.value}{  an approximated  p-value of the test.  }
  \item{method}{ the character string "Generalized Shapiro-Wilk test for multivariate normality".  }
  \item{data.name}{ a character string giving the name of the data set.  }
}
\references{ Villasenor-Alva, J.A. and Gonzalez-Estrada, E. (2009). A generalization of Shapiro-Wilk's test for multivariate normality. \emph{Communications in Statistics: Theory and Methods}, \bold{38} 11, 1870-1883. \url{http://dx.doi.org/10.1080/03610920802474465}}
\author{ Elizabeth Gonzalez-Estrada \email{egonzalez@colpos.mx}, Jose A. Villasenor-Alva}
\seealso{ \code{\link{shapiro.test}} and \code{\link{normal_test}}  for testing univariate normality.}
\examples{
# Example 1:  Testing multivariate normality on iris.virginica

# iris.virginica contains a set of measurements corresponding to 
# Iris virginica of famous  iris data set.

iris.virginica <- as.matrix(iris[iris$Species == "virginica", 1:4], ncol = 4) 
mvshapiro_test(iris.virginica)    


# Example 2:  Testing multivariate normality on the goats dataset
data(goats)
mvshapiro_test(as.matrix(goats))
}
\keyword{ htest }
\keyword{ multivariate analysis}