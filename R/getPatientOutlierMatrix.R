#
# Gavin a patMutMatrix, compute the corresponding patOutMatrix
#
getPatientOutlierMatrix <- function(patExpMatrix, th=2)
{
  expSd   <- apply(patExpMatrix, 2, sd)
  id <- expSd > 0
  expSd <- expSd[id]
  patExpMatrix <- patExpMatrix[, id]
  
  expMean <- apply(patExpMatrix, 2, mean)
	
  num <- dnorm(x=t(patExpMatrix), mean=expMean, sd=expSd, log=T)
  num <- t(num)
  
  numSd <- dnorm( x=t(expMean+th*expSd), mean=expMean, sd=expSd, log=T )
  
  y <- rep(numSd, each=dim(patExpMatrix)[1])
  y <- matrix(y, nrow=dim(patExpMatrix)[1], ncol=dim(patExpMatrix)[2])
	
  patOutMatrix <- num <= y
  
  id <- colSums(patOutMatrix) > 0
  patOutMatrix <- patOutMatrix[, id]
  
  return(patOutMatrix)
}