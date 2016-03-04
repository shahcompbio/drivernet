## Returns a matrix which each entry[g1, g2] represents the number of "events" between a gene g1 and another gene g2.
## There is an event when 
##   (1) gene "g1" is mutated in a patient x, 
##   (2) gene "g2" is a neighbor of "g1" in the influence graph, and
##   (3) the expression level "g2" is an outlier in the patient x.
.buildAggregateBipartiteGraph <- function(G, patMutMatrix, patOutMatrix) {
    bip <- matrix(0, nrow=dim(patMutMatrix)[2], ncol=dim(patOutMatrix)[2])
    rownames(bip) <- colnames(patMutMatrix)
    colnames(bip) <- colnames(patOutMatrix)
   
    drivers <- colnames(patMutMatrix)
    i <- 1
    for (i in 1:length(drivers)) {
        ## name of genes which are neighbor of "g" (drivers[i]) in the influence graph
        neighbor <- colnames(G)[which(G[drivers[i], ] == 1)]
        
        ## patients in which "g" is mutated
        pats <- which(patMutMatrix[, drivers[i]] == 1)

        if (length(pats) == 0) {
            next
        }
        
        j <- 1
        for (j in 1:length(pats)) {
            ## For each patient that "g" is mutated, first find all its outlier expression.
            ## pat_events is the intersect of these outliers and the neighbor of "g"
            pat_events <- intersect(neighbor, colnames(patOutMatrix)[which(patOutMatrix[pats[j], ] == 1)])
            if (length(pat_events) > 0) {
                bip[drivers[i], pat_events] <- bip[drivers[i], pat_events] + 1
            }
        }
    }
    
    ## post processing to remove disconnected (singelton) nodes
    count_mutated <- rowSums(bip)
    count_outlier <- colSums(bip)
    del_row <- which(count_mutated == 0)
    del_col <- which(count_outlier == 0)
    if (length(del_row) > 0) {
        bip <- bip[-del_row, ]
    }
    if (length(del_col) > 0) {
        bip <- bip[, -del_col]
    }    

    bip
}

## return a list of lists, which each of them gives the neighbors of the node
.neighborGraph <- function(G) {
  res <- list()
  rn <- rownames(G)
  cn <- colnames(G)
  
	## For each row, find all names of column that the corresponding entry is positive
  i <- 1
  for (i in 1:dim(G)[1]) {
      res[[rn[i]]] <- cn[which(G[i, ] > 0)]
  }
  res 
}

