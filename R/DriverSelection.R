.greedyGeneDriverSelection <- function(outf, patOutMatrix, patMutMatrix, influenceGraph, nG_out, fimatrix, maxNumOfDrivers) {
  ## If fimatrix is not provided, just use the binary patient mutation matrix
  ## If fimatrix is provided, we have to limit ourselves to use only patients & mutations existing in both matrices.
  if (!identical(fimatrix, NULL)) {
    pats <- intersect(rownames(patMutMatrix), rownames(fimatrix))
    cat(paste(dim(patMutMatrix)[1]-length(pats), "patients data removed because FI score not available\n"))
    patMutMatrix <- patMutMatrix[pats, ]
    patOutMatrix <- patOutMatrix[pats, ]
    fimatrix <- fimatrix[pats, ]
    patMutMatrix <- patMutMatrix[order(rownames(patMutMatrix)), ]
    patOutMatrix <- patOutMatrix[order(rownames(patOutMatrix)), ]
    fimatrix <- fimatrix[order(rownames(fimatrix)), ]
    
    muts <- intersect(colnames(patMutMatrix), colnames(fimatrix))
    cat(paste(dim(patMutMatrix)[2]-length(muts), "mutated genes removed because FI score not available\n"))
    
    patMutMatrix <- patMutMatrix[, muts]
    fimatrix <- fimatrix[, muts]
    patMutMatrix <- patMutMatrix[, order(colnames(patMutMatrix))]
    fimatrix <- fimatrix[, order(colnames(fimatrix))]
  }
  
  ## KEN: 2011-12-15: this total event is different from sum(patOutMatrix) (used below)... I do not have time to investigate why.
  totalEvents <- sum(.buildAggregateBipartiteGraph(influenceGraph, patMutMatrix, patOutMatrix))
  
  coveredEvents <- 0.0
  iter <- 0
  ## actualEvents keep track of in which patient, which outlier expression is explained by which mutation
  actualEvents <- list()
  res <- c()
  
  mutatedGenes <- colnames(patMutMatrix)
  mutationDegrees <- rep(0.0, dim(patMutMatrix)[2])
  names(mutationDegrees) <- mutatedGenes

  while (coveredEvents < totalEvents) {
    degRes <- .computeMutatedGeneDegrees(patMutMatrix, patOutMatrix, nG_out, fimatrix, mutatedGenes, mutationDegrees)
    mutationDegrees <- degRes$weightedDegrees
    actualDegrees <- degRes$actualDegrees
    maxMutDeg <- max(mutationDegrees)
    maxMutDegIndices <- which(mutationDegrees == maxMutDeg)
    actualMaxMutDeg <- actualDegrees[names(mutationDegrees)[maxMutDegIndices[[1]]]]

    ## Cannot cover any new event
    if (maxMutDeg == 0) {
      if (length(outf) > 0) {
        for (i in 1:length(outf)) {
          cat("exit prematurely!\n", file=outf[[i]], append=TRUE)
        }
      }
    break
    }

    maxMutDegGene <- names(mutationDegrees)[maxMutDegIndices[1]]
    if (length(outf) > 0) {
      for (i in 1:length(outf)) {
        cat("iter ", iter, " ) selected", maxMutDegGene, " ,  outliers covered by this driver: ", actualMaxMutDeg, "   total outliers explained:", sum(patOutMatrix), "\n", file=outf[[i]], append=T)
        if (length(maxMutDegIndices) > 1) {
            cat("   ties among : ", names(mutationDegrees)[maxMutDegIndices], "\n", file=outf[[i]], append=T)
        }
      }
    }        

    ## The max is removed
    coveredEvents <- coveredEvents + actualMaxMutDeg
    mutationDegrees <- mutationDegrees[-maxMutDegIndices[1]]
    updatedMatrices <- .deleteMutatedGene(patMutMatrix, patOutMatrix, nG_out[[maxMutDegGene]], maxMutDegGene)
    patMutMatrix <- updatedMatrices[[1]]
    patOutMatrix <- updatedMatrices[[2]]
    ## updatedMatrices[[3]] is the genesCovered list
    if (length(updatedMatrices[[3]]) > 1) {
      mutatedGenes <- colnames(patMutMatrix)[which( rowSums(influenceGraph[colnames(patMutMatrix), updatedMatrices[[3]]]) > 0)]
    } else {
      mutatedGenes <- colnames(patMutMatrix)[which( (influenceGraph[colnames(patMutMatrix), updatedMatrices[[3]]]) > 0)]
    }
    res <- c(res, maxMutDegGene)
    ## book-keeping for the outliers explained by this driver gene
    ## updatedMatrices[[4]] = those reset/cancelled outliers
    actualEvents[[maxMutDegGene]] <- updatedMatrices[[4]]
    rm(list=c("updatedMatrices"))
    gc() ##garbage collection
    iter <- iter + 1
    if (iter >= maxNumOfDrivers) {
      break
    }
  }
  list("drivers"=res, "actualEvents"=actualEvents, "totalEvents"=totalEvents)
}

## adjLists is the list of adjacency lists (usually called with nG_out, i.e. ordered by the rows in G)
## mutatedGenes: a list of mutation names
## mutation degree of these mutation
## returns a list of outlier event counts, value["g"] is the outlier event counts caused by mutation "g"
.computeMutatedGeneDegrees <- function(patMutMatrix, patOutMatrix, adjLists, fimatrix, mutatedGenes, weightedDeg) {
  weightedDeg[mutatedGenes] <- 0
  actualDeg <- weightedDeg
  
  i <- 1
  if (length(mutatedGenes) > 0) {
    for (i in 1:length(mutatedGenes)) {
      ## mutAdjList is the particular adjacency list of that gene 
      mutAdjList <- adjLists[[mutatedGenes[i]]]
      ## pats = patients in which "g" is mutated
      pats <- which(patMutMatrix[, mutatedGenes[i]] == 1)
      
      actualDeg[mutatedGenes[i]] <- sum(patOutMatrix[pats, intersect(mutAdjList, colnames(patOutMatrix))])
      
      if (identical(fimatrix, NULL)) {
        ## if fimatrix is not provided, the weightedDeg is the same as actual Deg
        ## This is a lot faster because of vectorization
        weightedDeg <- actualDeg
      } else {
        ## There should be a potential performance improvement here if we use vectorization.
        ## Accessing the matrix entries in this way is slow.
        if (length(pats) > 0) {
          for (j in 1:length(pats)) {
            ## weightedDeg["g"] = total number of outliers events of neighbor genes of "g" in all the patients whose gene "g" is mutated
            weightedDeg[mutatedGenes[i]] <- weightedDeg[mutatedGenes[i]] + fimatrix[pats[j], mutatedGenes[i]]*sum(patOutMatrix[pats[j], intersect(mutAdjList, colnames(patOutMatrix))])
          }
        }
      }
    }
  } else {
    weightedDeg[] <- 0
    actualDeg <- weightedDeg
  } 
  list(weightedDegrees=weightedDeg, actualDegrees=actualDeg)
}

## input: deleteMutatedGene(patOutMatrix, patMutMatrix, nG_out[[g]], g)
## g is the current choice
## remove that given mutation from patMutMatrix, all events covered by it are set to non-outliers
.deleteMutatedGene <- function(patMutMatrix, patOutMatrix, mutAdjList, geneName) {
    # mutatedGeneIndex = column index of geneName in patMutMatrix
    mutatedGeneIndex <- which(colnames(patMutMatrix) == geneName)
    # all patients that g is mutated
    pats <- which(patMutMatrix[, mutatedGeneIndex] == 1)
    j <- 1
    covered <- 0
    genesCovered <- c()
    actual_events_pat <- c()
    actual_events_gen <- c()
    ## for each patient...
    for (j in 1:length(pats)) {
        ## coveredEvents = set of outlier expression covered by the mutation of p of this patient
        coveredEvents <- intersect(mutAdjList, colnames(patOutMatrix)[which(patOutMatrix[pats[j], ] > 0)])
        if (length(coveredEvents) > 0) {
            ## set these coveredEvents to be not outlier
            patOutMatrix[pats[j], coveredEvents] <- 0
            ## actual_events_pat and _gen keep track of those "removed" outlier expression...
            actual_events_pat <- c(actual_events_pat, rep(rownames(patMutMatrix)[pats[j]], length(coveredEvents)))
            actual_events_gen <- c(actual_events_gen, coveredEvents)
            covered <- covered + length(coveredEvents)
            ## genesCovered keeps track of the covered events (not mentioning the patient though??)
            genesCovered <- c(genesCovered, coveredEvents)
        }
    }
    patMutMatrix <- patMutMatrix[, -mutatedGeneIndex]
    ## cbind: combine by column
    list(patMutMatrix, patOutMatrix, genesCovered, cbind(actual_events_pat, actual_events_gen))
}


