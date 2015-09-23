\name{goats}
\alias{goats}
\docType{data}
\title{
Zoometric measurements of creole goats
}
\description{
Zoometric measurements of 27 week old creole goats. These data were collected by Dorantes-Coronado (2013) as part of a preliminary study of the creole goat production system conducted in two locations in Mexico in order to start a genetic improvement program. Only one location is considered.

}
\usage{data(goats)}
\format{
  A data frame with 52 rows and 7 columns containing measurements (in kilograms and centimeters) on the following variables.
  \describe{
    \item{\code{body.weight}}{ }
    \item{\code{body.length}}{ }
    \item{\code{trunk.length}}{}
    \item{\code{withers.height}}{}
    \item{\code{thoracic.perimeter}}{}
    \item{\code{hip.length}}{}
    \item{\code{ear.length}}{}
  }
}

\source{
Dorantes-Coronado (2013).
}
\references{
Dorantes-Coronado, E.J. (2013). Estudio preliminar para el establecimiento de un programa de mejoramiento genetico de cabras en el Estado de Mexico. Tesis de doctorado. Colegio de Postgraduados, Mexico.
}
\examples{
data(goats)
plot(goats)
}
\keyword{datasets}
