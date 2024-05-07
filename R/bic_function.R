#' @title Calculate BIC values and generate the corresponding clustering results
#' @description \code{CalculateBIC} generates a matrix list of local clustering estimates and the corresponding BIC values based on the candidate lambda2 and lambda3 sequences.
#'
#' Dimensions of various arguments are as follows:
#' \itemize{
#' \item{n is the number of individuals}
#' \item{p = \code{nknots} + \code{order}, the number of basis functions}
#' }
#'
#' @param times A list with n elements, each of which contains a sequence of time points.
#' @param oridata A list with n elements, each of wich contains a sequence of observations corresponding to the times.
#' @param rangeval A numeric vector of length 2 defining the interval over which the functional data object can be evaluated; default value is c(0,1).
#' @param nknots An integer specifying the number of equally spaced knots; default value is 30.
#' @param order An integer specifying the order of B-splines, which is one higher than their degree. The default of 4 gives cubic splines.
#' @param nu A non-negative integer defining an order of a derivative. The default of 2 corresponds to measuring the roughness of a function by its integrated curvature.
#' @param tau A positive double parameter that controls the concavity of mcp function. The default value is 3.
#' @param K0 The pre-specified number of clusters; default value is 6.
#' @param rep_num The number of replicates for K-means; default value is 100.
#' @param lambda1 The tuning parameter for roughness of estimated functions.
#' @param lambda2_seq A sequence of tuning parameters for sparsity, that is, used to bring the estimated individual functions closer to the estimated center functions.
#' @param lambda3_seq A sequence of tuning parameters for clustering, that is, used to encourage similar estimated center functions to be the same, automatically identifying the number of clusters.
#' @param kappa A positive penalty parameter for ADMM algorithm. The default value is 1.
#' @param eps_outer The convergence tolerance and its default value is 0.0001.
#' @param max_iter The maximum number of iterations.
#'
#' @return A matrix, each element of which is a list of the following items.
#' \itemize{
#' \item{\code{funlc_result} A list of local clustering results.}
#' \item{\code{bic_val} The corresponding bic value.}
#' }
#'
#' @seealso \code{\link{LFC.fit}} for local clustering.
#'
#' @import fda
#' @import Matrix
#' @import Matrix
#' @import MASS
#' @import igraph
#'
#' @export
#'


CalculateBIC <- function(times, oridata, rangeval = c(0, 1), nknots = 30, order = 4, nu = 2, tau = 3, K0 = 6, rep_num = 100,
                         kappa = 1, eps_outer = 0.0001, max_iter = 100, lambda1, lambda2_seq, lambda3_seq){

  n <- length(oridata)
  p <- nknots + order
  L <- nknots + 1
  r <- order - 1
  N <- 0

  basisobj <- create.bspline.basis(rangeval = rangeval, nbasis = p, norder = order)

  Z_lst <- vector(mode = "list", length = n)
  for (i in 1:n) {
    Z_lst[[i]] <- eval.basis(times[[i]], basisobj)
    N <- N + nrow(Z_lst[[i]])
  }

  D <- eval.penalty(basisobj = basisobj, Lfdobj = nu)

  nlambda2 <- length(lambda2_seq)
  nlambda3 <- length(lambda3_seq)

  funlc.mat.lst <- matrix(list(), nrow = nlambda2, ncol = nlambda3)

  z <- 0
  total.num <- nlambda2 * nlambda3
  for (l1 in 1:nlambda2) {
    for (l2 in 1:nlambda3) {
      term1 <- rep(0, n)
      term2 <- rep(0, n)
      error <- rep(0, n)

      funlc_result <- LFC.fit(times = times, oridata = oridata, rangeval = rangeval, nknots = nknots, order = order,
                              nu = nu, tau = tau, K0 = K0, rep_num = rep_num, lambda1 = lambda1, lambda2 = lambda2_seq[l1],
                              lambda3 = lambda3_seq[l2], kappa = kappa, eps_outer = eps_outer, max_iter = max_iter)


      cls_no <- funlc_result$cls_num
      Beta <- funlc_result$Beta

      error_vec <- rep(0, n)
      df_vec <- rep(0, n)
      for (i in 1:n) {
        error_vec[i] <- sum((oridata[[i]] - as.numeric(Z_lst[[i]] %*% Beta[,i]))^2)
        ni <- length(oridata[[i]])
        Hi <- Z_lst[[i]] %*% solve(t(Z_lst[[i]]) %*% Z_lst[[i]] + lambda1*D) %*% t(Z_lst[[i]])
        df_vec[i] <- sum(diag(Hi))
      }
      error <- log(sum(error_vec)/N)
      df <- sum(df_vec) *  cls_no/n

      bn <- log(log(n*p))
      bic_val <- error + bn*log(N)*df/N

      funlc.mat.lst[[l1, l2]] <- list(funlc_result = funlc_result, bic_val = bic_val)


      z <- z + 1
      print(paste0('Finishing: ', 100 * round(z / total.num, 4), '%'))
      #print(c(l1, l2, cls_no, bic_val))
    }
  }

  return(funlc.mat.lst)
}
