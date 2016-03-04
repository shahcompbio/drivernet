\name{getPatientOutlierMatrix}
\alias{getPatientOutlierMatrix}
\title{
Compute the patient outlier matrix
}
\description{
Given a real-value patient expression matrix, compute the patient outlier matrix and remove the genes which are not expressed.
}
\usage{
getPatientOutlierMatrix(patExpMatrix, th=2)
}
\arguments{
  \item{patExpMatrix}{
Patient Expression Matrix
}
  \item{th}{
The threshold, the default value is 2 ( the gene expression values outside the 2 standard deviation range are considered as outliers).
}
}
\value{
  The output patOutMatrix is an input to the computeDrivers method 
}
\author{
Ali Bashashati, Reza Haffari, Jiarui Ding, Gavin Ha,  Kenneth Liu, Jamie Rosner and Sohrab Shah

Maintainer: Jiarui Ding <jiaruid@cs.ubc.ca>
}
\seealso{
\code{computeDrivers}
}
\examples{
data(samplePatientExpressionMatrix)
samplePatientOutlierMatrix = getPatientOutlierMatrix(samplePatientExpressionMatrix)
}
\keyword{ methods }
