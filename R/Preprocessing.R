## This hides all the options from the users because we want to keep the function easy to use.
preprocessMatrices <- function(patMutMatrix, patOutMatrix, influenceGraph) {
  .preprocessMatricesHelper(patMutMatrix, patOutMatrix, influenceGraph, reduceAll=FALSE, reducePatMut=TRUE, reducePatOut=TRUE, keepAllPatients=FALSE, reduceInfGraph=FALSE)
}

## by default reduceAll is TRUE and it overrides the default for all other options
## for finer options, set reduceAll to FALSE and choose the other options
.preprocessMatricesHelper <- function(patMutMatrix, patOutMatrix, influenceGraph, reduceAll = TRUE, reducePatMut = FALSE, reducePatOut = FALSE, keepAllPatients = TRUE, reduceInfGraph = FALSE) {
  if (reduceAll) {
      reduceInfGraph <- TRUE
      reducePatMut <- TRUE
      reducePatOut <- TRUE
      keepAllPatients <- FALSE
  }
  
  ## save intersection for reuse
  if (reducePatMut || reduceInfGraph)	mutInfluenceIntersect <- intersect(rownames(influenceGraph), colnames(patMutMatrix))
  if (reducePatOut || reduceInfGraph)	outInfluenceIntersect <- intersect(colnames(influenceGraph), colnames(patOutMatrix))
  
  ## patMutMatrix only keeps columns(mutation) that have corresponding rows in the influence graph
  if (reducePatMut) 
      patMutMatrix <- patMutMatrix[, mutInfluenceIntersect]
	
	## patOutMatrix only keeps columns(expression) that have corresponding columns in the influence graph																			
	if (reducePatOut)
		patOutMatrix <- patOutMatrix[, outInfluenceIntersect]

	## patientsIntersect = patients present in both mutation and expression data
	if (!keepAllPatients) {
		patientsIntersect <- intersect(rownames(patMutMatrix), rownames(patOutMatrix))	
		patMutMatrix <- patMutMatrix[patientsIntersect, ]
		patOutMatrix <- patOutMatrix[patientsIntersect, ]		
	}

	## influenceGraph only keeps the rows and columns of influenceGraph useful for this data set
	if (reduceInfGraph) influenceGraph <- influenceGraph[mutInfluenceIntersect, outInfluenceIntersect]
    
  list(patMutMatrix, patOutMatrix, influenceGraph)   
}

.calFunctionalImpactScore <- function(data, normalized=FALSE) {
  if (!is.data.frame(data))
  {
      stop("Type of parameter 'data' must be data frame.")
  }
  if (ncol(data) != 3)
  {
      stop("Parameter 'data' must have three columns.")
  }
  if (!is.numeric(data[, 3]))
  {
      stop("The third column of parameter 'data' must be numeric.")
  }
  patientID <- unique(data[,1])
  geneID <- unique(data[,2])
  result <- matrix(data=0.0, nrow=length(patientID), ncol=length(geneID))
  rownames(result) <- patientID
  colnames(result) <- geneID
  for (i in 1:nrow(data)) {
    if (!is.na(data[[i,3]])) {
      pID <- as.character(data[[i,1]])
      gID <- as.character(data[[i,2]])
      result[pID, gID] <- result[pID, gID] + data[[i,3]]
    }
  }
  if (normalized) {
    for (i in 1:nrow(result)) {
      result[i, ] <- result[i, ]/sum(result[i, ])
    }
  }
  result
} 

## apply logistic function to the FI matrix so that it can be used in the probabilistic version.
## let f(x) = 1 / (1 + e^(w*x + b))
## we want to find w and b such that 
## f(max(matrix)) = minBound
## f(min(matrix)) = maxBound
## then apply f to every element in the matrix
.applyLogisticFunction <- function(matrix, maxBound, minBound) {
  w <- log((1 / maxBound) - 1, base=exp(1)) - log((1 / minBound) - 1, base=exp(1))
  w <- w / (min(matrix) - max(matrix))
  b <- log((1 / maxBound) - 1, base=exp(1)) - min(matrix) * w
  resultMatrix <- 1/ (1 + exp(b + w * matrix))
  resultMatrix
}


## KEN: The following three functions are under development and potentially useless. Feel free to remove it in the future.
.makeAdjLists <- function(patMutMatrix, patOutMatrix, influenceGraph) {
  tmp <- .preprocessMatricesHelper(patMutMatrix, patOutMatrix, influenceGraph)
  pat_mut <- tmp[[1]]
  pat_out <- tmp[[2]]
  G <- tmp[[3]]
  rm(list=c("tmp"))
  gc() ##garbage collection
  
  neighLists <- .neighborGraph(G)
  
  adjLists <- vector("list", nrow(pat_mut))
  names(adjLists) <- rownames(pat_mut)
  for (i in 1:nrow(pat_mut)) {
    mutGenes <- colnames(pat_mut)[which(pat_mut[i, ]>0)]
    if (length(which(pat_mut[i, ]>0)) != 0) {
      adjLists[[i]] <- vector("list", length(which(pat_mut[i, ]>0)))
      names(adjLists[[i]]) <- mutGenes
      for (j in 1:length(which(pat_mut[i, ]>0))) {
        x <- intersect(neighLists[[mutGenes[j]]], colnames(pat_out)[which(pat_out[i, ]>0)])
        if (length(x) == 0) {
        } else {
          adjLists[[i]][[j]] <- x
        }      
      }
    }
  }
  
  perMutAdjList <- vector("list", ncol(pat_mut))
  names(perMutAdjList) <- colnames(pat_mut)
  for (i in 1:ncol(pat_mut)) {
    perMutAdjList[[i]] <- rownames(pat_mut)[which(pat_mut[, i]==1)]
  }
  
  list(al=adjLists, pal=perMutAdjList)
}

#.sumAdjList <- function(al) {
#  res = c()
#  for(i in 1:length(al$pal)) {
#    sum <- 0
#    if (length(al$pal[[i]])>0) {
#      for (j in 1:length(al$pal[[i]])) {
#        sum = sum + length(al$al[[al$pal[[i]][j]]][[colnames(pat_mut)[i]]])
#      }
#    }
#    res = c(res, sum)
#  }
#  res
#}

.newGreedy <- function(al) {
  
}