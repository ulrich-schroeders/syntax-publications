# title         Account for testlet structures
# GitHub        https://github.com/ulrich-schroeders/syntax-publications
# date          2016-11-06
# version       1.0.0
# reference     Schroeders, U., Robitzsch, A., & Schipolowski, S. (2014). A comparison of different psychometric approaches to modeling testlet structures: An example with c-Tests. Journal of Educational Measurement, 51(4), 400–418.

library(TAM)

# specify working directory, read in data
setwd("c:\\temp\\")
dat <- read.table( "text308.csv" , sep=";", dec=",", header=T , na="" )
testlet <- c(1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 4, 4, 4, 4, 4, 5, 5, 6, 6, 6)

# Rasch model
mod.rasch <- tam( resp = dat[,-1] , pid=dat[,1], 
	 control=list(maxiter=300, snodes=5000) )

# Testlet model
# testlets are dimensions, assign items to Q-matrix
TT <- length(unique(testlet))
Q <- matrix( 0 , nrow=ncol(dat[,-1]) , ncol= TT + 1)

Q[,1] <- 1 # First dimension constitutes g-factor
for (tt in 1:TT){ Q[ testlet == tt , tt+1 ] <- 1 }

# In a testlet model, all dimensions are uncorrelated among each other, 
# that is, all pairwise correlations are set to 0, which can be 
# accomplished with the "variance.fixed" command
library(combinat)
variance.fixed <- cbind( t(combn( TT+1,2 ) ) , 0 )

mod.testlet <- tam( resp = dat[,-1] , pid=dat[,1], Q = Q , 
	   variance.fixed = variance.fixed ,
	   control=list( snodes = 5000, QMC=FALSE, 
	   maxiter = 300 ) )

# Partial credit model
scores <- list()
testlet.names <- NULL
dat.pcm <- NULL
for (tt in 1:max(testlet)) {
	scores[[tt]] <- rowSums ( dat[,-1][, testlet == tt, drop=FALSE] )
	dat.pcm <- c(dat.pcm, list(c(scores[[tt]])))
	testlet.names <- append(testlet.names, paste("testlet",tt, sep=""))
}

dat.pcm <- as.data.frame(dat.pcm)
colnames(dat.pcm) <- testlet.names

mod.pcm <- tam( resp= dat.pcm , control=list( snodes = 5000 , QMC=FALSE , 	maxiter = 300 ) )

# Copula model
library(sirt)
mod.copula <- rasch.copula2 ( dat=dat[,-1] , itemcluster = testlet )
