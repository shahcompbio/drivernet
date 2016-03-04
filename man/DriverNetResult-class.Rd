\name{DriverNetResult-class}
\Rdversion{1.0.0}
\docType{class}
\alias{DriverNetResult-class}
\alias{actualEvents,DriverNetResult-method}
\alias{drivers,DriverNetResult-method}
\alias{totalEvents,DriverNetResult-method}

\title{Class \code{"DriverNetResult"}}
\description{
A class storing the result from running the \code{computeDrivers} function. It contains the list of driver mutations found, the detailed events covered by these drivers and the total number of events in the test data.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("DriverNetResult", drivers=..., 
actualEvents=..., totalEvents=...)}.
}
\section{Slots}{
  \describe{
    \item{\code{drivers}:}{Object of class \code{"character"} ~~ }
    \item{\code{actualEvents}:}{Object of class \code{"list"} ~~ }
    \item{\code{totalEvents}:}{Object of class \code{"numeric"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{actualEvents}{\code{signature(x = "DriverNetResult")}: ... }
    \item{drivers}{\code{signature(x = "DriverNetResult")}: ... }
    \item{totalEvents}{\code{signature(x = "DriverNetResult")}: ... }
	 }
}
\author{
Ali Bashashati, Reza Haffari, Jiarui Ding, Gavin Ha,  Kenneth Liu, Jamie Rosner and Sohrab Shah

Maintainer: Jiarui Ding <jiaruid@cs.ubc.ca>
}
\seealso{
\code{\link{computeDrivers}}
}
\examples{
data(sampleDriversList)
drivers(sampleDriversList)[1:10]
}
\keyword{classes}
