## Summarize the result returned by computeDrivers and computeRandomizedResult.
##
## If randResult is NULL, it means no randomized drivers are provided to calculate p-values.
##
## outputFolder = NULL -> no output any file
## outputFolder = "" -> current folder
resultSummary <- function (mainResult, randResult, patMutMatrix, influenceGraph, outputFolder = NULL, printToConsole = FALSE) {
  if (!inherits(mainResult, "DriverNetResult"))
  {
      stop("Type of parameter mainResult must be DriverNetResult.")
  }
  noRandResult <- identical(randResult, NULL)
  if (!noRandResult) {
    if ((!is.list(randResult)) || (!is.numeric(randResult[[1]])))
    {
        stop("Type of parameter randResult must be a list of numeric vectors.")
    }  
  }
  
  ## We want the different output (summary & gene_freq) to share the same timestamp
  timeString <- format(Sys.time(), "%Y_%m_%d_%H%M%S")
  summaryFileNames <- c()
  if (!identical(outputFolder, NULL)) {
    summaryFileNames <- c(paste(outputFolder, "summary_", timeString, ".txt", sep=""))
  }
  ## using "" as file name means printing to console
  if (printToConsole) {
    summaryFileNames <- c(summaryFileNames, "")
  }  
  
  actualEventsList <- actualEvents(mainResult)
  totalEventsCount <- totalEvents(mainResult)
  driversList <- drivers(mainResult)
  driversCount <- length(driversList)

  if (!noRandResult)
  {
    ## pValuesResult is a N x 1 matrix, rownames are driver names, N = numOfDriversFound(mainResult)    
    pValuesResult <- .computePValues(mainResult, randResult)
  }

  ## counts is a list of driver & count of mutated patients
  patientCounts <- colSums(patMutMatrix)[driversList]
  
  ## nodeDeg keeps the number of "influenced gene" of the drivers based on the influence graph
  ## order of drivers is the same as the list of drivers
  nodeDeg <- rowSums(influenceGraph[driversList, ])
  names(nodeDeg) <- NULL

	coveredEventsCount <- c()
  caseStrings <- c()
	casesCount <- c()

  for (i in 1:driversCount) {
    ## each element is the number of events covered by that driver
		coveredEventsCount <- c(coveredEventsCount, nrow(actualEventsList[[i]]))
    ## a list of patients (unique)
    uniqueCases <- unique(actualEventsList[[i]][, 1])
    ## a list of string, each string is the patients having this driver mutation
  	caseStrings <- c(caseStrings, paste(uniqueCases, collapse=","))
    ## contains the number of cases that have that have mutation in this driver gene
		casesCount <- c(casesCount, length(uniqueCases))
	}

  ## We only need to work on a smaller subset matrix of the influence graph, this speeds up the access to influence graph.
  ## This part of the code finds the drivers that are connected through the influence graph
  influenceGraph2 <- influenceGraph[driversList, driversList]
	connectedDrivers <- character(driversCount)
	for (u1 in 1:driversCount) {
    connectedNodes <- which(influenceGraph2[driversList[u1], ] == 1)
    ## The results are numbers (ranks of related drivers) instead of gene names. It is similar to neighborGraph() but it only deals with the drivers.
		connectedDrivers[u1] <- paste(connectedNodes, collapse=",")
	}
  
  ## If no random drivers given, the p-values are "N/A"
  if (noRandResult)
  {
    pValuesResult <- rep("N/A", driversCount)
  }

  ## combine the vectors by column to form the result summary table
  summaryTable <- cbind(seq(length(patientCounts)), driversList, patientCounts[driversList], pValuesResult, totalEventsCount, coveredEventsCount, nodeDeg, casesCount, connectedDrivers, caseStrings)
	colnames(summaryTable) <- c("rank", "gene", "mut","p-value", "total_events", "covered_events", "node_degree", "no.of.cases", "connected.drivers", "case")

  ## Write result summary to file (with name summaryFileNames)
  if (length(summaryFileNames) > 0) {
    for (i in 1:length(summaryFileNames)) {
      write.table(summaryTable, file=summaryFileNames[[i]], sep="\t", quote=FALSE, row.names=FALSE)
      if (!(summaryFileNames[[i]]=="")) {
        message(paste("Summary file written to: ", summaryFileNames[[i]], "\n", sep=""))
      }
    }
  }  

  if ((!noRandResult) & (!identical(outputFolder, NULL))) {
    write.table(summaryTable[pValuesResult<=0.01, c("gene", "no.of.cases")], file=paste(outputFolder, "gene_freq_", timeString, ".txt", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
    message(paste("Gene frequency file written to: ", outputFolder, "gene_freq_", timeString, ".txt", "\n", sep=""))
  }
  
  summaryTable
}

## Helper function to compute p-values from the randomized driver results 
.computePValues <- function(mainResult, randResult) {
  ## drivers found in the main run are now stored in main_sol
  main_sol <- drivers(mainResult)
  actualEventsList <- actualEvents(mainResult)
    
	counts <- c()
  for (i in 1:length(randResult)) {
    counts <- c(counts, randResult[[i]])
  }

  ## calculate p-values here
	pvals <- c()
	covered.net <- c() ## final length of this sould be length(main_sol)
	for (k in 1:length(main_sol)) {
		covered.net[k] <- nrow(actualEventsList[[k]])  ## nrow(actualEventsList[[k]]) is the number of events covered
    ## find where the covered.net sits with respect to counts
    ## find the faction of random drivers covering more events than this driver
    ## length(which(counts>covered.net[k])) -> number of drivers found in all runs which cover more events then the current driver
    ## length(counts) -> number of drivers found in all runs, thus the counts <- c(counts,0) above may make a difference
		pvals <- c(pvals, length(which(counts>covered.net[k])) / length(counts)) 
	}
  ## now pvals is a N x 1 matrix, rownames are driver names, N should be length(main_sol)
	pvals <- matrix(pvals, length(pvals), 1)
	rownames(pvals) <- main_sol
	return(pvals)
}