#' @title Create adjacency matrix from the pairwise difference matrix
#' @description \code{CreateAdjacency} creates an n-by-n sparse adjacency matrix from the matrix of centroid differences.
#'
#' @param V Matrix of centroid differences
#' @param n Number of points to cluster
#' @export



CreateAdjacency <- function(V, n) {
  differences <- apply(V,2,FUN=function(x) {norm(as.matrix(x),'f')})
  connected.ix <- which(differences == 0);
  index = t(combn(n,2));
  i <- index[connected.ix,1]
  j <- index[connected.ix,2]
  A <- Matrix(0, nrow = n, ncol = n, sparse = TRUE)
  A[(j-1)*n + i] <- 1
  return(A)
}

