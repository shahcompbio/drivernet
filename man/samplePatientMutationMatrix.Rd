\name{samplePatientMutationMatrix}
\alias{samplePatientMutationMatrix}
\docType{data}
\title{
Sample patient mutation matrix
}
\description{
A binary matrix representing gene mutations in patients.
}
\usage{data(samplePatientMutationMatrix)}
\format{
A two-dimensional binary matrix. The row names are patients, and collumn names are genes. If gene j is mutated in patient i, samplePatientMutationMatrix[i, j]=1.
}
\examples{
data(samplePatientMutationMatrix)
}
\keyword{datasets}
