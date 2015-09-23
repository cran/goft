\name{ig.fit}
\alias{ig.fit}
\title{ Fitting the Inverse Gaussian distribution to data}
\description{
Maximum likelihood estimators of the Inverse Gaussian distribution based on a random sample.
}
\usage{
ig.fit(x)
}
\arguments{
  \item{x}{ numeric data vector containing a random sample of positive real numbers.}
}

\value{
  Location and scale parameter estimates.
}
\author{ Elizabeth Gonzalez-Estrada \email{egonzalez@colpos.mx}, Jose A. Villasenor-Alva }
\seealso{ \code{\link{ig.test}} for testing the Inverse Gaussian distribution hypothesis.}
\examples{
x <-  rgamma(50,shape = 10)      
ig.fit(x)             # fitting an Inverse Gaussian distribution to x                        
}

