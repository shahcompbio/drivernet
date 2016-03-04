\name{DriverNet-package}
\alias{DriverNet-package}
\alias{DriverNet}
\docType{package}
\title{Drivernet: uncovering somatic driver mutations modulating transcriptional networks in cancer}
\description{
Description:  DriverNet is a package to predict functional important driver genes in cancer by integrating genome data (mutation and copy number variation data) and transcriptome data (gene expression data). The different kinds of data are combined by an influence graph, which is a gene-gene interaction network deduced from pathway data. A greedy algorithm is used to find the possible driver genes, which may mutated in a larger number of patients and these mutations will push the gene expression values of the connected genes to some extreme values. 
}
\details{
\tabular{ll}{
Package: \tab DriverNet\cr
Type: \tab Package\cr
Version: \tab 1.0.0\cr
Date: \tab 2012-03-21\cr
License: \tab GPL-3\cr
LazyLoad: \tab yes\cr
Depends: \tab methods\cr
}
Users would call preprocess_matrices on the three matrices to remove unnecessary data. The matrices returned are then be passed to other functions.
}
\author{
Ali Bashashati, Reza Haffari, Jiarui Ding, Gavin Ha,  Kenneth Liu, Jamie Rosner and Sohrab Shah

Maintainer: Jiarui Ding <jiaruid@cs.ubc.ca>
}
\references{
Drivernet: uncovering somatic driver mutations modulating transcriptional networks in cancer, in preparation. 
}
\keyword{ package }
\examples{

data(samplePatientMutationMatrix)
data(samplePatientOutlierMatrix)
data(sampleInfluenceGraph)
data(sampleGeneNames)

driversList = computeDrivers(samplePatientMutationMatrix, samplePatientOutlierMatrix, 
sampleInfluenceGraph, outputFolder=NULL, printToConsole=FALSE)

drivers(driversList)[1:10]

randomDriversResult = computeRandomizedResult(patMutMatrix=samplePatientMutationMatrix, 
patOutMatrix=samplePatientOutlierMatrix, influenceGraph=sampleInfluenceGraph, 
geneNameList= sampleGeneNames, outputFolder=NULL, printToConsole=FALSE, 
numberOfRandomTests=20, weight=FALSE, purturbGraph=FALSE, purturbData=TRUE)

res = resultSummary(driversList, randomDriversResult, samplePatientMutationMatrix, 
sampleInfluenceGraph, outputFolder=NULL, printToConsole=FALSE)}
