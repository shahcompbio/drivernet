## Usage ##
First open 
R, in Linux and Mac, just type R in shell. In R, type the following commands:
```R
library(DriverNet)

data(samplePatientMutationMatrix)
data(samplePatientOutlierMatrix)
data(sampleInfluenceGraph)
data(sampleGeneNames)

# The main function to compute drivers
driversList = computeDrivers(samplePatientMutationMatrix, samplePatientOutlierMatrix,
sampleInfluenceGraph, outputFolder=NULL, printToConsole=FALSE)

drivers(driversList)[1:10]

# random permute the gene labels to compute p-values
randomDriversResult = computeRandomizedResult(patMutMatrix=samplePatientMutationMatrix,
patOutMatrix=samplePatientOutlierMatrix, influenceGraph=sampleInfluenceGraph,
geneNameList= sampleGeneNames, outputFolder=NULL, printToConsole=FALSE,
numberOfRandomTests=20, weight=FALSE, purturbGraph=FALSE, purturbData=TRUE)

# Summarize the results
res = resultSummary(driversList, randomDriversResult, samplePatientMutationMatrix,
sampleInfluenceGraph, outputFolder=NULL, printToConsole=FALSE)
```

## Reproduce the paper results ##
First download the data file [paperData.tar.gz](ftp://ftp.bcgsc.ca/public/shahlab/DriverNet/paperData.tar.gz) and decompress it to the desired folder.  
Open R, and load the influence graph
```R
load("influenceGraph.rda") 
Then load one of the data file, e.g.,
load("GBM_data.rda")
```

and run

```R
driversList = computeDrivers(patMutMatrix, patOutMatrix,
influenceGraph, outputFolder=NULL, printToConsole=FALSE) 
```
will give you the rank list of genes.