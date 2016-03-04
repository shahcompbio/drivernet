\name{samplePatientOutlierMatrix}
\alias{samplePatientOutlierMatrix}
\docType{data}
\title{
Sample patient outlier matrix
}
\description{
A binary matrix representing gene expressions in patients.
}
\usage{data(samplePatientOutlierMatrix)}
\format{
  A two-dimensional binary matrix. The row names are patients, and collumn names are genes. If genes j is an outlier in patient i,  samplePatientOutlierMatrix[i, j] = 1.
}
\examples{
data(samplePatientOutlierMatrix)
}
\keyword{datasets}
