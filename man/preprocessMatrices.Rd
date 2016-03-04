\name{preprocessMatrices}
\alias{preprocessMatrices}
\title{
Remove unnecessary entries from matrices
}
\description{
Remove patients (rows) which are not in both matrices from patMutMatrix and patOutMatrix.
Remove mutations (columns) from patMutMatrix if they are not a row of influenceGraph.
Remove expressions (columns) from patOutMatrix if they are not a column of influenceGraph.
}
\usage{
preprocessMatrices(patMutMatrix, patOutMatrix, influenceGraph)
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
}
\value{
  \item{1 }{Reduced version of Patient Mutation Matrix}
  \item{2 }{Reduced version of Patient Expression(Outlier) Matrix}
  \item{3 }{Influence Graph Matrix}
}
\author{
Ali Bashashati, Reza Haffari, Jiarui Ding, Gavin Ha,  Kenneth Liu, Jamie Rosner and Sohrab Shah

Maintainer: Jiarui Ding <jiaruid@cs.ubc.ca>
}
\note{
The dimension of influenceGraph is not changed.
}
\examples{
data(samplePatientMutationMatrix)
data(samplePatientOutlierMatrix)
data(sampleInfluenceGraph)
tmp = preprocessMatrices(samplePatientMutationMatrix, samplePatientOutlierMatrix, 
sampleInfluenceGraph)

smallerPatientMutationMatrix = tmp[[1]]
smallerPatientOutlierMatrix = tmp[[2]]
influenceGraph = tmp[[3]]
}
\keyword{ methods }
