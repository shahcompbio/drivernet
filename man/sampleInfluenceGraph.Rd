\name{sampleInfluenceGraph}
\alias{sampleInfluenceGraph}
\docType{data}
\title{
Sample influence graph
}
\description{
A binary matrix representing influence graph between genes.
}
\usage{data(sampleInfluenceGraph)}
\format{
  A two-dimensional binary matrix. The row and column names are gene names. If two genes i and j are in the same pathway,  sampleInfluenceGraph[i, j] = 1.
}

\examples{
data(sampleInfluenceGraph)
}
\keyword{datasets}
