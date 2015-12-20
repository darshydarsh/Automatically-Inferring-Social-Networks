setwd(".../facebook-graph-data")

library(igraph)
library(expm)
library(sets)


# square root of a matrix
mat.sqrt <- function(A) {
  ei <- eigen(A)
  d <- ei$values
  d <- (d+abs(d))/2
  d2 <- sqrt(d)
  sqrt <- ei$vectors %*% diag(d2) %*% t(ei$vectors)
  return(sqrt)
}

# square-root-inverse of a matrix
mat.sqrt.inv <- function(A) {
  ei <- eigen(A)
  d <- ei$values
  d <- (d + abs(d))/2
  d2 <- 1/sqrt(d)
  d2[d == 0] <- 0
  sqrt.inv <- ei$vectors %*% diag(d2) %*% t(ei$vectors)
  return(sqrt.inv)
}

my.expm <- function(A, n) {
  A <- as.matrix(A)
  E <- diag(nrow(A))
  for (i in 1:n) {
    En <- (A %^% i)/factorial(i)
    E <- E + En
  }
  return(E)
}

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

spectral.clust <- function(M) {
  n <- dim(M)[1]
  A <- matrix(0, nrow=n, ncol=n)
  sigma <- 1
  for (i in 1:n) {
    for (j in 1:n) {
      A[i, j] <- exp(-sum((M[i, ] - M[j, ])^2)/(2*sigma^2))
    }
  }
  
  D <- matrix(0, nrow=n, ncol=n)
  for (i in 1:n) {
    D[i, i] <- sum(A[i, ])
  }
  
  sqrt.inv <- mat.sqrt.inv(D)
  L <- sqrt.inv %*% A %*% sqrt.inv
  
  eigen <- eigen(L)
  k <- k.largest.eigen(eigen$values, 0.7)
  X <- eigen$vectors[, 1:k]
  
  Y <- X
  for (i in 1:nrow(Y)) {
    sum <- sum(X[i, ]^2)
    if (sum == 0) {
      Y[i, ] <- X[i, ]
    } else {
      Y[i, ] <- X[i, ]/sum
    }
  }
  
  return(Y)
}

edit.dist <- function(x, y) {
  return(length(setdiff(union(x, y), intersect(x, y))))
}

Jaccard.ind <- function(x, y) {
  return(length(intersect(x, y))/length(union(x, y)))
}

Jaccard.dist <- function(x, y) {
  return(1 - Jaccard.ind(x, y))
}

my.sim <- function(x, y) {
  return(length(intersect(x, y))/min(length(x), length(y)))
}

cluster.density <- function(Adj, cl) {
  num.edges <- sum(Adj[cl, cl])/2
  total <- choose(length(cl), 2)
  density <- num.edges/total
  return(density)
}

discard.sparse <- function(clusters, rho) {
  density <- rep(0, length(clusters))
  for (i in 1:length(density)) {
    density[i] <- cluster.density(Adj, clusters[[i]])
  }
  
  l <- length(clusters)
  j <- 0
  for (i in 1:l) {
    j <- j + 1
    if (density[i] <= rho) {
      clusters[[j]] <- NULL
      j <- j - 1
    }
  }
  return(clusters)
}

recursive.remove <- function(clusters, rho) {
  for (i in 1:length(clusters)) {
    cat(paste("-------- Removal --------"))
    cat(paste('\n', "Cluster ", i, ":", sep=""))
    cat(paste('\n', "  Size before: ", old.len <- length(clusters[[i]]), sep=""))
    cat(paste('\n', "  Init density: ", cluster.density(Adj, clusters[[i]]), sep=""))
    while (cluster.density(Adj, clusters[[i]]) <= rho) {
      Adj.r <- Adj[clusters[[i]], clusters[[i]]]
      least.c <- which(rowSums(Adj.r) == min(rowSums(Adj.r)))
      clusters[[i]] <- clusters[[i]][-least.c]
    }
    cat(paste('\n', "  Size after: ", new.len <- length(clusters[[i]])))
    cat(paste('\n', "  removed: ", old.len - new.len, sep = ""))
    cat(paste('\n', "  Final density: ", cluster.density(Adj, clusters[[i]]), sep=""))
    cat(paste('\n'))
  }
  return(clusters)
}

forward.add <- function (clusters, rho) {
  for (i in 1:length(clusters)) {
    cat(paste("-------- Addition --------"))
    cat(paste('\n', "Cluster ", i, ":", sep=""))
    cat(paste('\n', "  Size before: ", old.len <- length(clusters[[i]]), sep=""))
    cat(paste('\n', "  Init density: ", cluster.density(Adj, clusters[[i]]), sep=""))
    while (cluster.density(Adj, clusters[[i]]) > rho) {
      Adj.r <- Adj
      Adj.r[clusters[[i]], ] <- 0
      Adj.r <- Adj.r[, clusters[[i]]]
      max.c <- which(rowSums(Adj.r) == max(rowSums(Adj.r)))
      if (max(rowSums(Adj.r)) > 0) {
        clusters[[i]] <- append(clusters[[i]], max.c)
        if (cluster.density(Adj, clusters[[i]]) < rho) {
          clusters[[i]] <- clusters[[i]][-max.c]
          break
        }
      } else {
        break
      }
    }
    cat(paste('\n', "  Size after: ", new.len <- length(clusters[[i]])))
    cat(paste('\n', "  Added: ", new.len - old.len, sep = ""))
    cat(paste('\n', "  Final density: ", cluster.density(Adj, clusters[[i]]), sep=""))
    cat(paste('\n'))
  }
  return(clusters)
}


# User 107 was returning an error at the get.adjacency() function:
# "Invalid (negative) vertex id", although I couldn't find any negative vertex
# id. I have removed 107 for now until I can figure out what's going on.
#
users <- c("0","348","414","686","698","1684","1912","3437","3980")
runtime <- rep(0, length(users))
Jaccard <- rep(0, length(users))
ed <- rep(0, length(users))
for (user in 1:length(users)) {
  start.time <- proc.time()
  circles.file <- file(paste(users[user], ".circles", sep=""))
  assign(paste(user, ".circles", sep=""), strsplit(readLines(circles.file), '\t'))
  close(circles.file)
  
  egofeat <- as.matrix(read.table(paste(users[user],".egofeat",sep=""), head=F))
  
  fr.feat <- as.matrix(read.table(paste(users[user],".feat",sep=""), head=F))
  fr.names <- fr.feat[, 1]
  fr.feat <- fr.feat[, -1]
  
  Adj <- read.table(paste(users[user],".edges",sep=""), head=F)
  Adj <- get.adjacency(graph.edgelist(as.matrix(Adj), directed=F))
  Adj <- Adj[fr.names, fr.names]
  Adj <- Adj/2
  
  E <- expm(as.matrix(Adj))
  YE <- spectral.clust(as.matrix(E))
  YF <- spectral.clust(fr.feat)
  
  if (length(fr.names) > 350) {
    k <- 10
  } else {
    k <- 6
  }
  
  cl.E <- kmeans(YE, centers=k)
  cl.F <- kmeans(YF, centers=k)
  
  clusters <- list(0)
  for (i in 1:k) {
    assign(paste(user, ".E.cluster", i, sep=""), which(cl.E$cluster==i))
    clusters [[i]] <- get(paste(user, ".E.cluster", i, sep=""))
  }
  
  for (i in 1:k) {
    assign(paste(user, ".F.cluster", i, sep=""), which(cl.F$cluster==i))
    clusters [[k + i]] <- get(paste("686.F.cluster", i, sep=""))
  }
  
    J.sim <- matrix(0, nrow=k, ncol=k)
  for (i in 1:k) {
    for (j in 1:k) {
      J.sim[i, j] <- Jaccard.ind(get(paste(user, ".E.cluster", i, sep="")), 
                                 get(paste(user, ".F.cluster", i, sep="")))
    }
  }
  Jaccard[user] <- sum(J.sim)/(k*k)

#   E.density <- rep(0, k)
#   F.density <- rep(0, k)
#   for (i in 1:k) {
#     E.density[i] <- cluster.density(Adj, get(paste("E.cluster", i, sep="")))
#     F.density[i] <- cluster.density(Adj, get(paste("F.cluster", i, sep="")))
#   }
  
#   clusters <- discard.sparse(clusters, 0.1)
#   clusters <- recursive.remove(clusters, 0.3)
#   clusters <- forward.add(clusters, 0.25)

  
#   cat(paste('\n', "User: ", users[user], "  # friends = ", length(fr.names), 
#             sep=""))
#   cat(paste('\n', "  E cluster sizes: ", sep=""))
#   cat(paste(cl.E$size))
#   cat(paste('\n', "      densities: ", sep=""))
#   cat(paste(round(E.density, 2)))
#   cat(paste('\n', "  F cluster sizes: ", sep=""))
#   cat(paste(cl.F$size))
#   cat(paste('\n', "      densities: ", sep=""))
#   cat(paste(round(F.density, 2)))
#   cat('\n')
#   print(round(J.sim, 2))
#   cat(paste('\n', "------------xxxxxxxx------------", sep=""))

  runtime[user] <- as.numeric(proc.time()[3] - start.time[3])


}

circles.file <- file("698.circles")
circles <- strsplit(readLines(circles.file), '\t')
close(circles.file)
# "circles" is a list of numerics, with the first element of each numeric the 
# circle lable e.g. "circle0" from the .circles file. The rest of the elements 
# of each numeric are the users in that circle.

egofeat <- as.matrix(read.table("698.egofeat", head=F))

fr.feat <- as.matrix(read.table("698.feat", head=F))
fr.names <- fr.feat[, 1]
fr.feat <- fr.feat[, -1]


Adj <- as.matrix(read.table("698.edges", head=F))
Adj <- get.adjacency(graph.edgelist(as.matrix(Adj), directed=F))
Adj <- Adj[fr.names, fr.names]
Adj <- Adj/2

E <- expm(as.matrix(Adj))

YE <- spectral.clust(as.matrix(E))
YF <- spectral.clust(fr.feat)

k <- 6
cl.E <- kmeans(YE, centers=k)
cl.F <- kmeans(YF, centers=k)



# Similarities between E clusters and F clusters
# cl.sim <- matrix(0, nrow=k, ncol=k)
# for (i in 1:k) {
#   for (j in 1:k) {
#     # Tried a lot of metrics here: Jaccard.ind, Jaccard.dist, edit.dist, 
#     # and size of intersection.
#     # 
#     # Inconclusive results as far as matching clusters goes.
#     # Some interesting results having to do with sizes of clusters.
#     # For user 0, E produced many small clusters and 1 very large cluster.
#     # One of the F clusters matched almost 96% with this very large cluster.
#     # Could be promising w.r.t. nesting.
#     cl.sim[i, j] <- Jaccard.ind(get(paste("E.cluster", i, sep="")), 
#                            get(paste("F.cluster", j, sep="")))
#   }
# }

clusters <- list(0)
for (i in 1:k) {
  assign(paste("686.", "E.cluster", i, sep=""), which(cl.E$cluster==i))
  clusters [[i]] <- get(paste("686.E.cluster", i, sep=""))
}

for (i in 1:k) {
  assign(paste("686.", "F.cluster", i, sep=""), which(cl.F$cluster==i))
  clusters [[k + i]] <- get(paste("686.F.cluster", i, sep=""))
}

clusters <- discard.sparse(clusters, 0.1)
clusters <- recursive.remove(clusters, 0.3)
clusters <- forward.add(clusters, 0.25)

for (user in 1:length(users)) {
  circles.file <- file(paste(users[user], ".circles", sep=""))
  assign(paste(user, ".circles", sep=""), strsplit(readLines(circles.file), '\t'))
  close(circles.file)
}

