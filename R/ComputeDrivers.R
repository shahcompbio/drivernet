## Comput a list of drivers ranked by the number of events covered by them
## This hides the use of fimatrix (functional impact matrix) from the users of the package because the implementation has not been finalized
## outputFolder = NULL -> no output any file
## outputFolder = "" -> current folder
computeDrivers <- function(patMutMatrix, patOutMatrix, influenceGraph, outputFolder = NULL, printToConsole = FALSE, weighted = FALSE) {
  .computeDriversHelper(patMutMatrix, patOutMatrix, influenceGraph, fimatrix = NULL, outputFolder, printToConsole, weighted)
}

## using fimatrix (Functional Impact Matrix)
.computeDriversHelper <- function(patMutMatrix, patOutMatrix, influenceGraph, fimatrix = NULL, outputFolder = NULL, printToConsole = FALSE, weighted = FALSE) {
  if (weighted) {
    stop("Weighted algorithm is not implemented.")
	}
  
  ## set maxNumOfDrivers to the number of mutations
  maxNumOfDrivers <- length(which(colSums(patMutMatrix)>0))
  outputFileNames <- c()
  if (!identical(outputFolder, NULL)) {
    outputFileNames <- c(paste(outputFolder, "compute_drivers", format(Sys.time(), "%Y_%m_%d_%H%M%S"), ".txt", sep=""))
  }
  if (printToConsole) {
    outputFileNames <- c(outputFileNames, "")
  }
  
  ## nG_out = genes whose expression is affected by dna mutation in g
  nG_out <- .neighborGraph(influenceGraph[intersect(colnames(patMutMatrix), rownames(influenceGraph)), intersect(colnames(influenceGraph), colnames(patOutMatrix))])
  gc() ##garbage collection

  tmp <- .preprocessMatricesHelper(patMutMatrix, patOutMatrix, influenceGraph)
	patMutMatrix <- tmp[[1]] 
	patOutMatrix <- tmp[[2]]
	influenceGraph <- tmp[[3]] 
  rm(list=c("tmp"))
  gc() ##garbage collection

	if (weighted == FALSE) {
    res <- .greedyGeneDriverSelection(outputFileNames, patOutMatrix, patMutMatrix, influenceGraph, nG_out, fimatrix, maxNumOfDrivers)
	}

  if (!identical(outputFolder, NULL)) {
    message(paste("Log file written to: ", outputFileNames[[1]], "\n", sep=""))
  }
  
  ## Package the result as a DriverNetResult object
  driverNetRes <- new ("DriverNetResult", drivers=res$drivers, actualEvents=res$actualEvents, totalEvents=res$totalEvents)
}

