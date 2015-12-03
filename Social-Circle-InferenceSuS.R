setwd(".../facebook-graph-data")

library(igraph)
library(Matrix)
library(expm)

# Square Root of a Matrix
mat.sqrt <- function(A) {
  ei <- eigen(A)
  d <- ei$values
  d <- (d+abs(d))/2
  d2 <- sqrt(d)
  sqrt <- ei$vectors %*% diag(d2) %*% t(ei$vectors)
  return(sqrt)
}

# Square-root-inverse Of a Matrix
mat.sqrt.inv <- function(A) {
  ei <- eigen(A)
  d <- ei$values
  d <- (d + abs(d))/2
  d2 <- 1/sqrt(d)
  d2[d == 0] <- 0
  sqrt.inv <- ei$vectors %*% diag(d2) %*% t(ei$vectors)
  return(sqrt.inv)
}

# obtaining the Largest Eigenvalue
k.largest.eigen <- function(values, threshold) {
  k <- 0
  sum <- sum(values^2)
  var <- 0
  for (i in 1:length(values)) {
    if (var < threshold*sum) {
      var <- var + values[i]^2
      k <- k + 1
    } else {
      return(k)
    }
  }
  return(k)
}

#Applying Spectral Clustering
spectral.clust <- function(M) {
  n <- dim(M)[1]
  A <- matrix(0, nrow=n, ncol=n)
  sigma <- 1
  for (i in 1:n) {
    for (j in 1:n) {
      A[i, j] <- exp(-sum((M[i, ] - M[j, ])^2)/(2*sigma^2))
    }
  }
  
  #Matrix with sums of each row
  D <- matrix(0, nrow=n, ncol=n)
  for (i in 1:n) {
    D[i, i] <- sum(A[i, ])
  }
  
  #Square root of inverse Matrix
  sqrt.inv <- mat.sqrt.inv(D)
  L <- sqrt.inv %*% A %*% sqrt.inv
  
  eigen <- eigen(L)
  k <- k.largest.eigen(eigen$values, 0.7)
  X <- eigen$vectors[, 1:k]
  
#Summing over each row
  Y <- X
  for (i in 1:nrow(Y)) {
    Y[i, ] <- X[i, ]/sum(X[i, ]^2)
  }
  
  return(Y)
}

#Adjacency Matrix
Adj <- read.table("0.edges", head=F)
Adj <- get.adjacency(graph.edgelist(as.matrix(Adj), directed=F))
Adj <- Adj/2

#Exponential Adjacency Matrix
E <- expm(as.matrix(Adj))

#Reading in Files
egofeat <- as.matrix(read.table("0.egofeat", head=F))

fr.feat <- as.matrix(read.table("0.feat", head=F))
fr.names <- fr.feat[, 1]
fr.feat <- fr.feat[, -1]


#Applying Spectral Clustering
Y <- spectral.clust(as.matrix(E))

Y1 <- spectral.clust(fr.feat)
