\name{o3}
\alias{o3}
\docType{data}
\title{
Ozone concentrations over 0.165 ppm
}
\description{
Ozone levels in parts per million (ppm) registered in the southwest  of Mexico City from January 1, 2008 to April 30th, 2016, which exceed the 0.165  ppm threshold.
}
\usage{data(o3)}
\format{
    A data frame with 47 observations on the following 2 variables.
  \describe{
    \item{\code{date}}{date in the format yyyy-mm-dd}
    \item{\code{ozone_level}}{a numeric vector containing the ozone levels above 0.165 ppm}
  }
}

\source{
Atmospheric Monitoring Office of Mexico City.
}
\examples{
data(o3)
plot(o3)
}
\keyword{datasets}
