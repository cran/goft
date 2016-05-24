\name{o3max}
\alias{o3max}
\docType{data}
\title{
Monthly maximum ozone concentrations
}
\description{
Maximum ozone concentrations (in ppm) per month registered in the southwest  of Mexico City from 2008 to 2015 during the months of March to May.   
}
\usage{data(o3max)}
\format{
  A data frame with 24 observations on the following 2 variables.
  \describe{
    \item{\code{month}}{months in which observations were registered.}
    \item{\code{o3max_level}}{a numeric vector containing the maximum ozone leveles per month.}
  }
}
\source{
Atmospheric Monitoring Office of Mexico City.  \url{http://www.aire.df.gob.mx/default.php?opc='aqBjnmU='}
}

\examples{
data(o3max)
plot(o3max)
}
\keyword{datasets}
