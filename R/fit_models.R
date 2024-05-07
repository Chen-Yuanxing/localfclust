#' @title Local clustering for functional data
#' @description \code{LFC.fit} performs local clustering on functional data via proximal average ADMM algorithm.
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
#' @param lambda2 The tuning parameter for sparsity, that is, used to bring the estimated individual functions closer to the estimated center functions.
#' @param lambda3 The tuning parameter for clustering, that is, used to encourage similar estimated center functions to be the same, automatically identifying the number of clusters.
#' @param kappa A positive penalty parameter for ADMM algorithm. The default value is 1.
#' @param eps_outer The convergence tolerance and its default value is 0.0001.
#' @param max_iter The maximum number of iterations.
#'
#' @return A list of the following items.
#' \itemize{
#' \item \code{Beta}: A p * n matrix of B-spline coefficients for all individuals.
#' \item \code{Beta_ini}: A p * n matrix of initial smoothed B-spline coefficients for all individuals.
#' \item \code{centers}: A matrix of cluster centers.
#' \item \code{basisobj}: A B-spline basis object.
#' \item \code{cls_mem}: A vector of integers indicating the cluster to which each individual is allocated.
#' \item \code{cls_num}: The number of final clusters.
#' }
#'
#' @import fda
#' @import Matrix
#' @import MASS
#' @import igraph
#'
#' @export
#'
#'

LFC.fit <- function(times, oridata, rangeval = c(0, 1), nknots = 30, order = 4, nu = 2, tau = 3, K0 = 6, rep_num = 100,
                    lambda1, lambda2, lambda3, kappa = 1, eps_outer = 0.0001, max_iter = 100){

  ########################################################################
  ##################      The intermediate functions.     ################
  ########################################################################

  ### Function to update Beta ###
  update_Beta <- function(ZTY_lst, ZTZ_lst, Beta, Alpha, D, tau, varsigma1, lambda1, lambda2, iter_outer, max_iter_inner){

    p <- nrow(Beta)
    n <- ncol(Beta)

    U_middle <- matrix(0, nrow = p, ncol = n)
    V_middle <- matrix(0, nrow = p, ncol = n)
    Beta_old <- matrix(0, nrow = p, ncol = n)
    Beta_new <- matrix(0, nrow = p, ncol = n)

    for (i in 1:n) {

      eta_old <- 1
      U_middle[,i] <- Beta[,i]
      Beta_old[,i] <- U_middle[,i]

      iter_inner <- 1
      res_inner  <- 1
      eps_inner  <- 0.01/iter_outer

      while (res_inner > eps_inner & iter_inner <= max_iter_inner) {

        V_middle[,i]     <- U_middle[,i] - varsigma1 * as.numeric((ZTZ_lst[[i]] + lambda1*D) %*% U_middle[,i] - ZTY_lst[[i]])
        center_index_opt <- assign_center_label(V_middle[,i], Alpha)
        Beta_new[,i]     <- ProximalMCP(x = V_middle[,i] - Alpha[,center_index_opt], varsigma1, tau, lambda2) + Alpha[,center_index_opt]
        eta_new          <- (1 + sqrt(1 + 4*eta_old^2))/2
        U_middle[,i]     <- Beta_new[,i] + (eta_old - 1)/eta_new * (Beta_new[,i] - Beta_old[,i])


        res_inner    <- sqrt(sum((Beta_new[,i] - Beta_old[,i])^2))
        eta_old      <- eta_new
        iter_inner   <- iter_inner + 1
        Beta_old[,i] <- Beta_new[,i]
      }
      #print(iter_inner)
    }
    return(Beta_new)
  }

  ### Function to update Alpha ###
  update_Alpha <- function(Beta, Alpha, Psi, Gamma, Omega, tau, kappa, varsigma2, lambda2, iter_outer, max_iter_inner){

    p  <- nrow(Beta)
    n  <- ncol(Beta)
    K0 <- ncol(Alpha)

    Alpha_old <- Alpha
    Alpha_new <- Alpha
    u_middle  <- as.vector(Alpha)

    eta_old    <- 1
    iter_inner <- 1
    res_inner  <- 1
    eps_inner  <- 0.01/iter_outer


    OmegaTpsi   <- t(Omega) %*% as.vector(Psi)
    OmegaTgamma <- t(Omega) %*% as.vector(Gamma)/kappa
    OmegaTOmega <- t(Omega) %*% Omega

    while (res_inner > eps_inner & iter_inner <= max_iter_inner) {

      v_middle  <- u_middle - varsigma2*kappa/n * as.numeric(OmegaTOmega %*% u_middle - OmegaTpsi + OmegaTgamma)
      V_middle  <- matrix(v_middle, nrow = p, ncol = K0)
      index_mat <- matrix(0, nrow = n, ncol = K0)

      for (i in 1:n) {
        dis <- rep(0, K0)
        for (k in 1:K0) {
          dis[k] <- sqrt(sum((Beta[,i] - V_middle[, k])^2))
        }
        index_mat[i,which.min(dis)] <- 1
      }

      for (k in 1:K0) {
        if (sum(index_mat[,k] == 1) == 0) Alpha_new[,k] <- V_middle[,k]
        else {
          Alpha_new[,k] <- rep(0, p)
          for (i in 1:n) {
            if (index_mat[i,k] == 1) Alpha_new[,k] <- Alpha_new[,k] + ProximalMCP(x = V_middle[,k] - Beta[,i], varsigma2, tau, lambda2) + Beta[,i]
          }
          Alpha_new[,k] <- Alpha_new[,k]/sum(index_mat[,k] == 1)
        }
      }

      eta_new  <- (1 + sqrt(1 + 4*eta_old^2))/2
      u_middle <- as.vector(Alpha_new) + (eta_old - 1)/eta_new * as.vector(Alpha_new - Alpha_old)

      res_inner  <- norm(Alpha_new - Alpha_old, type = "F")
      eta_old    <- eta_new
      iter_inner <- iter_inner + 1
      Alpha_old  <- Alpha_new
    }


    return(Alpha_new)
  }

  ### Function to update Psi ###
  update_Psi <- function(Psi, Alpha, Gamma, index, n, L, r, tau, varsigma3, kappa, lambda3, iter_outer, max_iter_inner){

    p <- nrow(Psi)
    S <- nrow(index)

    U_middle <- matrix(0, nrow = p, ncol = S)
    V_middle <- matrix(0, nrow = p, ncol = S)
    Psi_old  <- matrix(0, nrow = p, ncol = S)
    Psi_new  <- matrix(0, nrow = p, ncol = S)

    for (s in 1:S) {
      eta_old <- 1
      U_middle[,s] <- Psi[,s]
      Psi_old[,s]  <- U_middle[,s]

      iter_inner <- 1
      res_inner  <- 1
      eps_inner  <- 0.01/iter_outer

      while (res_inner > eps_inner & iter_inner <= max_iter_inner) {

        V_middle[,s] <- U_middle[,s] - varsigma3*kappa*(U_middle[,s] - (Alpha[,index[s,1]] - Alpha[,index[s,2]] + Gamma[,s]/kappa)) / (n*L)
        Psi_new[,s]  <- rep(0, p)
        for (l in 1:L) {
          Psi_new[,s] <- Psi_new[,s] + ProximalAverageMCP(x = V_middle[,s], l, r, varsigma3, tau, lambda3)
        }
        Psi_new[,s]  <- Psi_new[,s]/L

        for (l in 1:L) {
          Psi_new[l:(l+r),s] <- ProximalMCP(Psi_new[l:(l+r),s], varsigma3, tau, lambda3)
        }

        eta_new      <- (1 + sqrt(1 + 4*eta_old^2))/2
        U_middle[,s] <- Psi_new[,s] + (eta_old - 1)/eta_new * (Psi_new[,s] - Psi_old[,s])
        res_inner    <- sqrt(sum((Psi_new[,s] - Psi_old[,s])^2))
        Psi_old[,s]  <- Psi_new[,s]
        eta_old      <- eta_new
        iter_inner   <- iter_inner + 1

      }
      #print(iter_inner)
    }
    return(Psi_new)
  }

  ### Function to update Gamma ###
  update_Gamma <- function(Gamma_old, Alpha, Psi, index, kappa){

    p <- nrow(Gamma_old)
    S <- ncol(Gamma_old)

    Gamma <- matrix(0, nrow = p, ncol = S)
    for (s in 1:S) {
      i <- index[s,1]
      j <- index[s,2]
      Gamma[,s] <- Gamma_old[,s] + kappa*(Alpha[,i] - Alpha[,j] - Psi[,s])
    }
    return(Gamma)
  }

  ### Function to assign center label to each individual  ###
  assign_center_label <- function(u, Alpha){

    K0 <- ncol(Alpha)

    dis_vec <- rep(0, K0)
    for (k in 1:K0) {
      dis_vec[k] <- sum((u - Alpha[,k])^2)
    }
    return(which.min(dis_vec))

  }

  ### Function to calculate varsigma1 ###
  cal_varsigma1 <- function(Z_lst, D, lambda){

    n <- length(Z_lst)
    Lip_vec <- rep(0, n)

    for (i in 1:n) {
      term <- t(Z_lst[[i]]) %*% Z_lst[[i]]
      eigen_result <- eigen(term)
      Lip_vec[i] <- max(as.numeric(eigen_result$values))
    }
    Lip <- max(Lip_vec) + lambda*max(as.numeric(eigen(D)$values))
    varsigma1 <- 1/Lip

    return(varsigma1)
  }

  ### Function to calculate varsigma2 ###
  cal_varsigma2 <- function(Omega, kappa) {
    Lip <- kappa*max(as.numeric(eigen(t(Omega) %*% Omega)$values))
    varsigma2 <- 1/Lip

    return(varsigma2)
  }


  #############################################################################################
  ##################      The main body of Proximal average ADMM algorithm     ################
  #############################################################################################

  index <- t(combn(K0, 2))
  S     <- nrow(index)
  n     <- length(oridata)
  p     <- nknots + order
  L     <- nknots + 1
  r     <- order - 1

  basisobj <- create.bspline.basis(rangeval = rangeval, nbasis = p, norder = order)

  Z_lst <- vector(mode = "list", length = n)
  Y_lst <- vector(mode = "list", length = n)
  for (i in 1:n) {
    Z_lst[[i]] <- eval.basis(times[[i]], basisobj)
    Y_lst[[i]] <- oridata[[i]]
  }
  D <- eval.penalty(basisobj = basisobj, Lfdobj = nu)

  coef_ini <- GenerateInitialCoefficients(times = times, oridata = oridata, rangeval = rangeval, nknots = nknots,
                                          order = order, nu = nu, lambda1 = lambda1, K0 = K0, rep_num = rep_num)

  Beta_old  <- coef_ini$Beta_ini
  Alpha_old <- coef_ini$Alpha_ini

  ZTY_lst <- vector(mode = "list", length = n)
  ZTZ_lst <- vector(mode = "list", length = n)
  for (i in 1:n) {
    ZTY_lst[[i]] <- t(Z_lst[[i]]) %*% Y_lst[[i]]
    ZTZ_lst[[i]] <- t(Z_lst[[i]]) %*% Z_lst[[i]]
  }


  Psi_old   <- matrix(0, nrow = p, ncol = S)
  Omega_lst <- vector(mode = "list", length = S)
  for (s in 1:S) {
    i <- index[s,1]
    j <- index[s,2]
    Psi_old[,s] <- Alpha_old[,i] - Alpha_old[,j]
    e <- rep(0, K0)
    e[i] <- 1
    e[j] <- -1
    Omega_lst[[s]] <- kronecker(t(e), diag(1,p))
  }
  Omega <- do.call("rbind", Omega_lst)

  Gamma_old <- matrix(0, nrow = p, ncol = S)


  varsigma1 <- cal_varsigma1(Z_lst, D, lambda1)
  varsigma2 <- cal_varsigma2(Omega, kappa)
  varsigma3 <- 1/kappa
  #print(c(varsigma1, varsigma2, varsigma3))
  if(tau <= varsigma1 | tau <= varsigma2 | tau <= varsigma3) print("error stepsize")

  res_outer <- 1
  iter      <- 1

  while (res_outer > eps_outer & iter <= max_iter) {
    Beta_new  <- update_Beta(ZTY_lst, ZTZ_lst, Beta_old, Alpha_old, D, tau, varsigma1, lambda1, lambda2, iter, max_iter_inner = 5)
    Alpha_new <- update_Alpha(Beta_new, Alpha_old, Psi_old, Gamma_old, Omega, tau, kappa, varsigma2, lambda2, iter, max_iter_inner = 5)
    Psi_new   <- update_Psi(Psi_old, Alpha_new, Gamma_old, index, n, L, r, tau, varsigma3, kappa, lambda3, iter, max_iter_inner = 5)
    Gamma_new <- update_Gamma(Gamma_old, Alpha_new, Psi_new, index, kappa)

    res_outer_vec <- rep(0, S)
    for (s in 1:S) {
      i <- index[s,1]
      j <- index[s,2]
      res_outer_vec <- sum((Alpha_new[,i] - Alpha_new[,j] - Psi_new[,s])^2)
    }
    res_outer <- sqrt(sum(res_outer_vec))
    res_outer <- norm(Beta_new - Beta_old, type = "F")/norm(Beta_old, type = "F")

    for (l in 1:L) {
      Ad_final <- CreateAdjacency(Psi_new[l:(l+r),], K0)

      if (sum(as.matrix(Ad_final)) != 0) {
        G_final <- graph.adjacency(Ad_final, mode = 'upper')
        cls_final <- components(G_final)

        for (g in 1:cls_final$no) {
          Alpha_new[l:(l+r),cls_final$membership == g] <- apply(as.matrix(Alpha_new[l:(l+r),cls_final$membership == g]), MARGIN = 1, mean)
        }
      }
    }

    Beta_old  <- Beta_new
    Alpha_old <- Alpha_new
    Psi_old   <- Psi_new
    Gamma_old <- Gamma_new
    iter      <- iter + 1
  }


  for (i in 1:n) {
    center_index_opt <- assign_center_label(Beta_new[,i], Alpha_new)
    Beta_new[,i] <- ProximalMCP(Beta_new[,i] - Alpha_new[,center_index_opt], varsigma1, tau, lambda2) + Alpha_new[,center_index_opt]
  }


  group <- rep(0, n)
  for (i in 1:n) {
    dis_vec <- rep(0, K0)
    for (k in 1:K0) {
      dis_vec[k] <- sum((Beta_new[,i] - Alpha_new[,k])^2)
    }
    group[i] <- which.min(dis_vec)
  }

  Beta_trans <- matrix(0, nrow = p, ncol = n)
  for (g in unique(group)) {
    Beta_trans[,group == g] <- Alpha_new[,g]
  }

  for (j in 1:p) {
    dist_mat <- as.matrix(dist(Beta_trans[j,]))
    adj_vec <- dist_mat[t(upper.tri(dist_mat))]
    adj_vec[abs(adj_vec) <= 0.2] <- 0
    adj_mat <- matrix(adj_vec, nrow = 1)
    Ad_final <- CreateAdjacency(adj_mat, n);
    G_final <- graph.adjacency(Ad_final, mode = 'upper')
    cls_final <- components(G_final)
    for (k in 1:cls_final$no) {
      Beta_trans[j, cls_final$membership == k] <- mean(Beta_trans[j,cls_final$membership == k])
    }
  }

  cls_num     <- length(unique(group))
  Alpha_trans <- matrix(0, nrow = p, ncol = cls_num)
  for (k in 1:cls_num) {
    Alpha_trans[,k] <- apply(as.matrix(Beta_trans[,group == unique(group)[k]]), MARGIN = 1, mean)
  }

  cls_mem <- rep(0, n)
  for (i in 1:n) {
    dis_vec <- rep(0, cls_num)
    for (k in 1:cls_num) {
      dis_vec[k] <- sum((Beta_trans[,i] - Alpha_trans[,k])^2)
    }
    cls_mem[i] <- which.min(dis_vec)
  }

  return(list(Beta = Beta_trans, Beta_ini = coef_ini$Beta_ini, centers = Alpha_trans, basisobj = basisobj,
              cls_mem = cls_mem, cls_num = cls_num))
}


