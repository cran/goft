\name{sn_test}
\alias{sn_test}
\title{ Test for skew normality }
\description{
Transformation test of fit for skew-normality  (Gonzalez-Estrada and Cosmes, 2019).
}
\usage{
sn_test(x, method = "transf")
}
\arguments{
  \item{x}{ a numeric data vector containing a random sample of  size n. }
  \item{method}{ a character string giving the name of the test to be used. So far the only option is \verb{"transf"}. }
}
\details{
The random sample is transformed to approximately normal variables. Shapiro-Wilk test  is used for testing normality.

The sample size should be larger than 50 in order to have a reasonable control of the type I error probability.
}
\value{
  A list with class \code{"htest"} containing the following components.
  \item{p.value}{  an approximated p-value of the test.  }
  \item{method}{ the character string "Shapiro-Wilk test for skew-normal distributions".  }
  \item{data.name}{ a character string giving the name of the data set.  }
}
\references{ Gonzalez-Estrada, E. and Cosmes, W. (2020). 
Shapiro-Wilk test for skew normal distributions based on data transformations.  \emph{Journal of Statistical Computation and Simulation}, \bold{89} 17, 3258-3272.
\url{https://doi.org/10.1080/00949655.2019.1658763}
}
\author{ Elizabeth Gonzalez-Estrada (\email{egonzalez@colpos.mx})}
\examples{
data(strength)                # loading the "strength" data set
y  <- strength$strain   
sn_test(y)                    # testing skew normality on the strain variable  
}
\keyword{ htest }



