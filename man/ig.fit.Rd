\name{ig.fit}
\alias{ig.fit}
\title{ Fitting the Inverse Gaussian distribution to data}
\description{
Maximum likelihood estimators for the Inverse Gaussian distribution based on a random sample.
}
\usage{
ig.fit(x)
}
\arguments{
  \item{x}{ a numeric data vector containing a random sample of positive real numbers.}
}

\value{
  Parameter estimates.
}
\author{ Elizabeth Gonzalez-Estrada \email{egonzalez@colpos.mx}, Jose A. Villasenor-Alva }
\seealso{ \code{\link{ig.test}} for testing the Inverse Gaussian distribution hypothesis.}
\examples{
x <-  rgamma(50, 10)   # simulating a random sample from the gamma distribution 
ig.fit(x)              # fitting an inverse Gaussian distribution to x                        
}

