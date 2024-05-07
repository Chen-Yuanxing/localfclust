#' @title Calculate GCV values from a series of candidate lambda1
#' @description \code{CalculateGCV} selects the optimal lambda1 by minimizing GCV.
#'
#'
#'
#' @param times A list with n elements, each of which contains a sequence of time points.
#' @param oridata A list with n elements, each of wich contains a sequence of observations corresponding to the times.
#' @param rangeval A numeric vector of length 2 defining the interval over which the functional data object can be evaluated; default value is c(0,1).
#' @param nknots An integer specifying the number of equally spaced knots; default value is 30.
#' @param order An integer specifying the order of B-splines, which is one higher than their degree. The default of 4 gives cubic splines.
#' @param nu A non-negative integer defining an order of a derivative. The default of 2 corresponds to measuring the roughness of a function by its integrated curvature.
#' @param lambda1_seq A sequence of candidate tuning parameters for roughness of estimated functions.
#'
#' @export
#'

CalculateGCV <- function(times, oridata, rangeval = c(0, 1), nknots = 30, order = 4, nu = 2, lambda1_seq){

  n <- length(oridata)
  p <- nknots + order

  basisobj <- create.bspline.basis(rangeval = rangeval, nbasis = p, norder = order)

  Z_lst <- vector(mode = "list", length = n)
  for (i in 1:n) {
    Z_lst[[i]] <- eval.basis(times[[i]], basisobj)
  }

  D <- eval.penalty(basisobj = basisobj, Lfdobj = nu)

  nlambda1 <- length(lambda1_seq)

  gcv_value <- rep(0, nlambda1)

  for (l in 1:nlambda1) {
    term1 <- rep(0, n)
    term2 <- rep(0, n)
    error <- rep(0, n)

    for (i in 1:n) {
      ni <- length(oridata[[i]])
      Hi <- Z_lst[[i]] %*% solve(t(Z_lst[[i]]) %*% Z_lst[[i]] + lambda1_seq[l]*D) %*% t(Z_lst[[i]])

      term1[i] <- mean((oridata[[i]] - as.numeric(Hi %*% oridata[[i]]))^2)
      term2[i] <- mean(diag(diag(1, ni) - Hi))^2
      error[i] <- term1[i]/term2[i]
    }

    gcv_value[l] <- sum(error)

  }
  lambda1_op <- lambda1_seq[which.min(gcv_value)]
  return(lambda1_op)

}
