\name{sampleDriversList}
\alias{sampleDriversList}
\docType{data}
\title{
Sample DriverNet result
}
\description{
An object of DriverNetResult class that obtained from running the code{computeDrivers} function.
}
\usage{data(sampleDriversList)}
\format{
Formal class 'DriverNetResult' with 3 slots (1) @ drivers: a character vector of the driver gene name; (2) @ actualEvents: the covered events in each mutated patient; (3) @totalEvents: the total number of outliers
}
\seealso{
\code{computeDrivers}
}
\examples{
data(sampleDriversList)

drivers(sampleDriversList)[1:10]
}
\keyword{datasets}
