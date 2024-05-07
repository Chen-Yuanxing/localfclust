#' @title Generate the initial smoothed B-spline coefficients and the centered coefficients
#' @description \code{GenerateInitialCoefficients} generates the initial B-spline coefficients.
#'
#' Dimensions of various arguments are as follows:
#' \itemize{
#' \item{n is the number of individuals}
#' \item{p = \code{nknots} + \code{order}, the number of basis functions}
#' }
#'
#'
#' @param times A list with n elements, each of which contains a sequence of time points.
#' @param oridata A list with n elements, each of wich contains a sequence of observations corresponding to the times.
#' @param rangeval A numeric vector of length 2 defining the interval over which the functional data object can be evaluated; default value is c(0,1).
#' @param nknots An integer specifying the number of equally spaced knots; default value is 30.
#' @param order An integer specifying the order of B-splines, which is one higher than their degree. The default of 4 gives cubic splines.
#' @param nu A non-negative integer defining an order of a derivative. The default of 2 corresponds to measuring the roughness of a function by its integrated curvature.
#' @param lambda1 The tuning parameter for roughness of estimated functions.
#' @param K0 The pre-specified number of clusters; default value is 6.
#' @param rep_num The number of replicates for K-means; default value is 100.
#'
#' @return \code{Beta_ini} The p * n matrix, whose columns correspond to individuals.
#' @return \code{Alpha_ini} The p * \code{K0} matrix, whose columns correspond to centers.
#'
#' @export
#'


GenerateInitialCoefficients <- function(times, oridata, rangeval = c(0, 1), nknots = 30, order = 4, nu = 2, lambda1 = 1e-5, K0 = 6, rep_num = 100){

  n <- length(oridata)
  p <- nknots + order

  basisobj <- create.bspline.basis(rangeval = rangeval, nbasis = p, norder = order)

  Z_lst <- vector(mode = "list", length = n)
  for (i in 1:n) {
    Z_lst[[i]] <- eval.basis(times[[i]], basisobj)
  }

  D <- eval.penalty(basisobj = basisobj, Lfdobj = nu)

  Beta_ini <- matrix(0, nrow = p, ncol = n)
  for (i in 1:n) {
    Beta_ini[,i] <- solve(t(Z_lst[[i]]) %*% Z_lst[[i]] + lambda1*D) %*% t(Z_lst[[i]]) %*% oridata[[i]]
  }


  error <- rep(0, rep_num)
  center_lst <- vector(mode = "list", length = rep_num)
  cluster_lst <- vector(mode = "list", length = rep_num)
  for (l in 1:rep_num) {
    ini_cluster <- kmeans(t(Beta_ini), centers = K0)
    center_lst[[l]] <- t(ini_cluster$centers)
    cluster_lst[[l]] <- ini_cluster$cluster

    obj_value <- rep(0, n)
    for (i in 1:n) {
      obj_value[i] <- sum((as.numeric(Z_lst[[i]] %*% center_lst[[l]][,cluster_lst[[l]][i]]) - oridata[[i]])^2)
    }
    error[l] <- mean(obj_value)
  }
  Alpha_ini <- center_lst[[which.min(error)]]

  return(list(Beta_ini = Beta_ini, Alpha_ini = Alpha_ini))
}



