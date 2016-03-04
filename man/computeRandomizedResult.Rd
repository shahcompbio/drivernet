\name{computeRandomizedResult}
\alias{computeRandomizedResult}
\title{
Randomly compute lists of driver mutations
}
\description{
It renames the mutations in patMutMatrix with a randomized list of gene names and computes a list of driver mutations.
It will repeat this process by numberOfRandomTests times.
The result could be passed to the resultSummary method to calculate p-value for the driver mutations.
}
\usage{
computeRandomizedResult(patMutMatrix, patOutMatrix, influenceGraph, 
geneNameList, outputFolder = NULL, printToConsole = FALSE, 
numberOfRandomTests = 500, weighted = FALSE, purturbGraph = FALSE, 
purturbData = TRUE)
}
\arguments{
  \item{patMutMatrix}{
Patient Mutation Matrix
}
  \item{patOutMatrix}{
Patient Expression(Outlier) Matrix
}
  \item{influenceGraph}{
Influence Graph Matrix
}
  \item{geneNameList}{
A list of gene names that the new mutation names will be randomly chosen from.
}
  \item{outputFolder}{
The folder to store the log. If set to NULL, no log files will be written. If set to "", the log will be written to the current folder.
}
  \item{printToConsole}{
If set to TRUE, progress and result of the function will be printed to the console.
}
  \item{numberOfRandomTests}{
Number of lists of randomized driver mutations to be generated.
}
  \item{weighted}{
Must be set to FALSE in this version.
}
  \item{purturbGraph}{
Must be set to FALSE in this version.
}
  \item{purturbData}{
Must be set to TRUE in this version.
}
}
\value{
A list of numeric vectors where each vector stores the number of events covered by the randomly identified drivers. The list can be passed to the resultSummary method to compute p-values.
}
\author{
Ali Bashashati, Reza Haffari, Jiarui Ding, Gavin Ha,  Kenneth Liu, Jamie Rosner and Sohrab Shah

Maintainer: Jiarui Ding <jiaruid@cs.ubc.ca>
}
\examples{
data(samplePatientMutationMatrix)
data(samplePatientOutlierMatrix)
data(sampleInfluenceGraph)
data(sampleGeneNames)

randomDriversResult = computeRandomizedResult(patMutMatrix=samplePatientMutationMatrix, 
patOutMatrix=samplePatientOutlierMatrix, influenceGraph=sampleInfluenceGraph, 
geneNameList=sampleGeneNames, outputFolder=NULL, printToConsole=FALSE, 
numberOfRandomTests=20, weight=FALSE, purturbGraph=FALSE, purturbData=TRUE)

%load(system.file("data", "sampleData.rda", package="DriverNet"))
%rand_result = computeRandomizedResult(reducedPatientMutationMatrix, reducedPatientOutlierMatrix, reducedInfluenceGraph, geneNames, outputFolder=NULL, numberOfRandomTests=5)
}
\keyword{ methods }
