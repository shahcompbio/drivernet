\name{computeDrivers}
\alias{computeDrivers}
\title{
Compute a list of driver mutations.
}
\description{
Use a greedy algorithm to rank a list of driver mutations.
}
\usage{
computeDrivers(patMutMatrix, patOutMatrix, influenceGraph, outputFolder = NULL, 
printToConsole = FALSE, weighted = FALSE)
}
\arguments{
  \item{patMutMatrix}{
Patient Mutation Matrix
}
  \item{patOutMatrix}{
Patient Outlier Matrix
}
  \item{influenceGraph}{
Influence Graph Matrix
}
  \item{outputFolder}{
The folder to store the log. If set to NULL, no log files will be written. If set to "", the log will be written to the current folder.
}
  \item{printToConsole}{
If set to TRUE, progress and result of the function will be printed to the console.
}
  \item{weighted}{
Must be set to FALSE in this version.
}
}
\value{
An object of DriverNetResult class that can be passed to the \code{resultSummary} method.
}
\author{
Ali Bashashati, Reza Haffari, Jiarui Ding, Gavin Ha,  Kenneth Liu, Jamie Rosner and Sohrab Shah

Maintainer: Jiarui Ding <jiaruid@cs.ubc.ca>
}
\examples{
data(samplePatientMutationMatrix)
data(samplePatientOutlierMatrix)
data(sampleInfluenceGraph)
driversList = computeDrivers(samplePatientMutationMatrix, samplePatientOutlierMatrix, 
sampleInfluenceGraph, outputFolder=NULL, printToConsole=FALSE)

drivers(driversList)[1:10]
}
\keyword{ methods }