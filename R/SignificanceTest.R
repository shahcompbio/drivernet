## Runs the greedy algorithm on randomized input to produce random result list of drivers for significance test
## outputFolder = NULL -> no output any file
## outputFolder = "" -> current folder
computeRandomizedResult <- function(patMutMatrix, patOutMatrix, influenceGraph, geneNameList, outputFolder = NULL, printToConsole = FALSE, numberOfRandomTests = 500,  weighted = FALSE, purturbGraph = FALSE, purturbData = TRUE) {
  if (weighted) {
    stop("Weighted algorithm is not implemented.")
  }
  if (purturbGraph) {
    stop("Purturb graph option is not implemented.")
  }
  if (purturbData == FALSE) {
    stop("Purturb data option must be TRUE at this moment.")
  }  
  ## set maxNumOfDrivers to the number of mutations
  maxNumOfDrivers <- length(which(colSums(patMutMatrix)>0))
  out_fname <- c()
  if (!identical(outputFolder, NULL)) {
    out_fname <- c(paste(outputFolder, "randomized_result_", format(Sys.time(), "%Y_%m_%d_%H%M%S"), ".txt", sep=""))	
  }
  if (printToConsole) {
    out_fname <- c(out_fname, "")
  }  
  
  if (purturbGraph == FALSE) {
      nG_out <- .neighborGraph(influenceGraph[intersect(colnames(patMutMatrix), rownames(influenceGraph)), intersect(colnames(influenceGraph), colnames(patOutMatrix))])   ## genes whose expression is affected by dna mutation in g
  }

  coverageResults <- vector(mode="list", length=numberOfRandomTests)
  i <- 1
  for (i in 1:numberOfRandomTests) {
      randomPatMutMatrix <- patMutMatrix
      randomPatOutMatrix <- patOutMatrix
      if (purturbData) {
          ## purturb the data by randomizing the gene names
          randomizedOutlierNames <- geneNameList[sample(1:length(geneNameList))[1:ncol(patOutMatrix)]] 
          randomizedMutationNames <- geneNameList[sample(1:length(geneNameList))[1:ncol(patMutMatrix)]]
          colnames(randomPatOutMatrix) <- randomizedOutlierNames
          colnames(randomPatMutMatrix) <- randomizedMutationNames
      }
      if (purturbGraph) {
        numColInfluenceGraph=ncol(influenceGraph) 
        ## purturb influenceGraph as well
        infGraphNewGeneNames <- geneNameList[sample(1:length(geneNameList))[1:numColInfluenceGraph]]
        colnames(influenceGraph) <- infGraphNewGeneNames
        rownames(influenceGraph) <- infGraphNewGeneNames
        nG_out <- .neighborGraph(influenceGraph[intersect(colnames(randomPatMutMatrix), rownames(influenceGraph)), intersect(colnames(influenceGraph), colnames(randomPatOutMatrix))])   ## genes whose expression is affected by dna mutation in g
      }
      ## reuse these intersections
      affectedGenesIntersection <- intersect(colnames(influenceGraph), colnames(randomPatOutMatrix))
      patientIntersection <- intersect(rownames(randomPatMutMatrix), rownames(randomPatOutMatrix))
      mutatedGenesIntersection <- intersect(rownames(influenceGraph), colnames(randomPatMutMatrix))
      ## pre-process the new matrices
      influenceGraph2 <- influenceGraph[mutatedGenesIntersection, affectedGenesIntersection]
      randomPatOutMatrix2 <- randomPatOutMatrix[patientIntersection, affectedGenesIntersection]
      randomPatMutMatrix2 <- randomPatMutMatrix[patientIntersection, mutatedGenesIntersection] 
      if (weighted) {
        ## Not implemented
        # drivers[[i]] <- .greedyGeneDriverSelection_weighted(out_fname, randomPatOutMatrix2, randomPatMutMatrix2, influenceGraph2, maxNumOfDrivers)[[1]]   
      } else {
        runResult <- .greedyGeneDriverSelection(out_fname, randomPatOutMatrix2, randomPatMutMatrix2, influenceGraph2, nG_out, NULL, maxNumOfDrivers)
        if (!is.null(runResult$drivers)) {
            len <- length(runResult$drivers)
            coverageVector <- vector(mode="integer", length=len)
            names(coverageVector) <- runResult$drivers
            if (len>0) {
                for (k in 1:len) {
                  ## nrow(actualEvents[[i]][[k]]) equals the number of events covered by the k-th driver found in the i-th run
                  coverageVector[[k]] <- nrow(runResult$actualEvents[[k]])
                }
            }
            coverageResults[[i]] <- coverageVector
        }
      }
  }
  
  if (!identical(outputFolder, NULL)) {
    message(paste("Log file written to: ", out_fname[[1]], "\n", sep=""))
  }  

  coverageResults
}