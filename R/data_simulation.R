#' @title Simulate data for local functional clustering
#' @description Generate synthetic data for Examples 1-4 as in the simulation study (Chen et al., 2024).
#'
#' @param seed An integer specifying the random seed used for random data generation.
#' @param setpoint1 A positive double parameter that sets the location of the startpoint for heterogeneous subinterval in Example 1.
#' @param n An integer specifying the number of functions.
#' @param n0 An integer specifying the number of repeated measurements.
#' @param sigma A positive double parameter that controls the variance of random errors


#' @return A list of the following items.
#' \itemize{
#' \item \code{location}: A list of n vector and each vector collects the observed time points of the sequence data.
#' \item \code{data.lst}: A list of n vectors and each vector collects the observed replications.
#' \item \code{clus.true}: True cluster membership vector.
#' }
#'
#' @import MASS
#' @export SimulateExample1
#' @export SimulateExample2
#' @export SimulateExample3
#' @export SimulateExample4
#'
#' @name SimulateData
NULL

#' @rdname SimulateData
SimulateExample1 <- function(seed, setpoint1 = 0.6, n = 100, n0 = 100, sigma = 0.75){

  set.seed(seed)

  Times <- seq(0, 1, len = n0)
  class1 <- 1:(n/2)
  class2 <- (n/2+1):n
  setpoint2 <- 1


  mu <- rep(0, n0)
  rho1 <- 0.3
  sig <- matrix(0, nrow = n0, ncol = n0)
  for (i1 in 1:n0) {
    for (i2 in 1:n0) {
      sig[i1, i2] <- sigma^2 * rho1^abs(i1-i2)
    }
  }

  times <- vector(mode = "list", length = n)
  oridata_list <- vector(mode = "list", length = n)

  for (i in 1:n) {
    times[[i]] <- Times

    if (i %in% class1) {
      c1 <- 0*times[[i]]
      for (l in 1:n0) {
        if (times[[i]][l] >= setpoint1 & times[[i]][l] <= setpoint2) c1[l] <- 1 + 1*sin(2*pi*(times[[i]][l] - setpoint1)/(setpoint2 - setpoint1))
        else c1[l] <- 1
      }
      epsilon <- mvrnorm(1, mu = mu, Sigma = sig)
      oridata_list[[i]] <- c1 + epsilon
    }


    if (i %in% class2) {
      c2 <- 0*times[[i]]
      for (l in 1:n0) {
        if (times[[i]][l] >= setpoint1 & times[[i]][l] <= setpoint2) c2[l] <- 1 - 1*sin(2*pi*(times[[i]][l] - setpoint1)/(setpoint2 - setpoint1))
        else c2[l] <- 1
      }

      epsilon <- mvrnorm(1, mu = mu, Sigma = sig)
      oridata_list[[i]] <- c2 + epsilon
    }
  }

  return(list(location = times, data_lst = oridata_list, clus.true = c(class1, class2)))
}

#' @rdname SimulateData
SimulateExample2 <- function(seed, n = 100, n0 = 100, sigma = 0.5){

  set.seed(seed)

  Times <- seq(0, 1, len = n0)
  class1 <- 1:(2*n/5)
  class2 <- (2*n/5+1):(7*n/10)
  class3 <- (7*n/10+1):n


  setpoint1 <- 0.2
  setpoint2 <- 0.6
  setpoint3 <- 1


  times <- vector(mode = "list", length = n)
  oridata_list <- vector(mode = "list", length = n)

  mu <- rep(0, n0)
  rho1 <- 0.3
  sig <- matrix(0, nrow = n0, ncol = n0)
  for (i1 in 1:n0) {
    for (i2 in 1:n0) {
      sig[i1, i2] <- sigma^2 * rho1^abs(i1-i2)
    }
  }



  for (i in 1:n) {
    times[[i]] <- Times

    if (i %in% class1) {
      c1 <- 0*times[[i]]
      for (l in 1:n0) {
        if (times[[i]][l] >= 0 & times[[i]][l] <= setpoint1) c1[l] <- 1 + 1.5*sin(1*pi*(times[[i]][l] - 0)/(setpoint1 - 0))
        else if (times[[i]][l] <= setpoint2) c1[l] <- 1
        else c1[l] <- 1 - 1.5*sin(1*pi*(times[[i]][l] - setpoint2)/(setpoint3 - setpoint2))
      }
      epsilon <- mvrnorm(1, mu = mu, Sigma = sig)
      oridata_list[[i]] <- c1 + epsilon
    }


    if (i %in% class2) {
      c2 <- 0*times[[i]]
      for (l in 1:n0) {
        if (times[[i]][l] >= 0 & times[[i]][l] <= setpoint1) c2[l] <- 1 + 1.5*sin(1*pi*(times[[i]][l] - 0)/(setpoint1 - 0))
        else if (times[[i]][l] <= setpoint2) c2[l] <- 1 - 1.5*sin(1*pi*(times[[i]][l] - setpoint1)/(setpoint2 - setpoint1))
        else c2[l] <- 1
      }
      epsilon <- mvrnorm(1, mu = mu, Sigma = sig)
      oridata_list[[i]] <- c2 + epsilon
    }


    if (i %in% class3) {
      c3 <- 0*times[[i]]
      for (l in 1:n0) {
        if (times[[i]][l] >= 0 & times[[i]][l] <= setpoint1) c3[l] <- 1 + 1.5*sin(1*pi*(times[[i]][l] - 0)/(setpoint1 - 0))
        else if (times[[i]][l] <= setpoint2) c3[l] <- 1 - 1.5*sin(1*pi*(times[[i]][l] - setpoint1)/(setpoint2 - setpoint1))
        else c3[l] <- 1 + 1.5*sin(1*pi*(times[[i]][l] - setpoint2)/(setpoint3 - setpoint2))
      }

      epsilon <- mvrnorm(1, mu = mu, Sigma = sig)
      oridata_list[[i]] <- c3 + epsilon
    }
  }


  return(list(location = times, data_lst = oridata_list, clus.true = c(class1, class2, class3)))
}

#' @rdname SimulateData
SimulateExample3 <- function(seed, n = 100, n0 = 100, sigma = 0.5){

  set.seed(seed)

  Times <- seq(0, 1, len = n0)
  class1 <- 1:(n/4)
  class2 <- (n/4+1):(n/2)
  class3 <- (n/2+1):(3*n/4)
  class4 <- (3*n/4+1):n


  setpoint1 <- 0.3
  setpoint2 <- 0.6
  setpoint3 <- 1


  times <- vector(mode = "list", length = n)
  oridata_list <- vector(mode = "list", length = n)

  mu <- rep(0, n0)
  rho1 <- 0.3
  sig <- matrix(0, nrow = n0, ncol = n0)
  for (i1 in 1:n0) {
    for (i2 in 1:n0) {
      sig[i1, i2] <- sigma^2 * rho1^abs(i1-i2)
    }
  }


  for (i in 1:n) {
    times[[i]] <- Times

    if (i %in% class1) {
      c1 <- 0*times[[i]]
      for (l in 1:n0) {
        if (times[[i]][l] <= setpoint1) c1[l] <- 1 + 1.5*sin(1*pi*(times[[i]][l] - 0)/(setpoint1 - 0))
        else if (times[[i]][l] <= setpoint2) c1[l] <- 1 - 1.5*sin(1*pi*(times[[i]][l] - setpoint1)/(setpoint2 - setpoint1))
        else c1[l] <- 1 + 1.5*sin(1*pi*(times[[i]][l] - setpoint2)/(setpoint3 - setpoint2))
      }
      epsilon <- mvrnorm(1, mu = mu, Sigma = sig)
      oridata_list[[i]] <- c1 + epsilon
    }

    if (i %in% class2) {
      c2 <- 0*times[[i]]
      for (l in 1:n0) {
        if (times[[i]][l] <= setpoint2) c2[l] <- 1 + 1.5*sin(1*pi*(times[[i]][l] - 0)/(setpoint1 - 0))
        else if (times[[i]][l] > setpoint2) c2[l] <- 1
        else c2[l] <- 1 - 1.5*sin(1*pi*(times[[i]][l] - setpoint2)/(setpoint3 - setpoint2))
      }
      epsilon <- mvrnorm(1, mu = mu, Sigma = sig)
      oridata_list[[i]] <- c2 + epsilon
    }

    if (i %in% class3) {
      c3 <- 0*times[[i]]
      for (l in 1:n0) {
        if (times[[i]][l] <= setpoint1) c3[l] <- 1
        else if (times[[i]][l] <= setpoint2) c3[l] <- 1 - 1.5*sin(1*pi*(times[[i]][l] - setpoint1)/(setpoint2 - setpoint1))
        else c3[l] <- 1
      }
      epsilon <- mvrnorm(1, mu = mu, Sigma = sig)
      oridata_list[[i]] <- c3 + epsilon
    }


    if (i %in% class4) {
      c4 <- 0*times[[i]]
      for (l in 1:n0) {
        if (times[[i]][l] <= setpoint1) c4[l] <- 1
        else if (times[[i]][l] <= setpoint2) c4[l] <- 1
        else c4[l] <- 1 + 1.5*sin(1*pi*(times[[i]][l] - setpoint2)/(setpoint3 - setpoint2))
      }

      epsilon <- mvrnorm(1, mu = mu, Sigma = sig)
      oridata_list[[i]] <- c4 + epsilon
    }
  }

  return(list(location = times, data_lst = oridata_list, clus.true = c(class1, class2, class3, class4)))
}

#' @rdname SimulateData
SimulateExample4 <- function(seed, n = 100, n0 = 100, sigma = 0.75){

  set.seed(seed)

  Times <- seq(0, 1, len = n0)
  class1 <- 1:(n/4)
  class2 <- (n/4+1):(n/2)
  class3 <- (n/2+1):(3*n/4)
  class4 <- (3*n/4+1):n


  setpoint1 <- 0.3
  setpoint1 <- 0.2
  setpoint2 <- 0.4
  setpoint3 <- 0.6
  setpoint4 <- 0.8
  setpoint5 <- 1

  times <- vector(mode = "list", length = n)
  oridata_list <- vector(mode = "list", length = n)

  mu <- rep(0, n0)
  rho1 <- 0.3
  sig <- matrix(0, nrow = n0, ncol = n0)
  for (i1 in 1:n0) {
    for (i2 in 1:n0) {
      sig[i1, i2] <- sigma^2 * rho1^abs(i1-i2)
    }
  }


  for (i in 1:n) {
    times[[i]] <- Times

    if (i %in% class1) {
      c1 <- 0*times[[i]]
      for (l in 1:n0) {
        if (times[[i]][l] <= setpoint1) c1[l] <- 1
        else if (times[[i]][l] <= setpoint2) c1[l] <- 1 + 1*sin(1*pi*(times[[i]][l] - setpoint1)/(setpoint2 - setpoint1))
        else if (times[[i]][l] <= setpoint3) c1[l] <- 1 - 1*sin(1*pi*(times[[i]][l] - setpoint2)/(setpoint3 - setpoint2))
        else if (times[[i]][l] <= setpoint4) c1[l] <- 1 + 1*sin(1*pi*(times[[i]][l] - setpoint3)/(setpoint4 - setpoint3))
        else c1[l] <- c1[l] <- 1
      }
      epsilon <- mvrnorm(1, mu = mu, Sigma = sig)
      oridata_list[[i]] <- c1 + epsilon
    }

    if (i %in% class2) {
      c2 <- 0*times[[i]]
      for (l in 1:n0) {
        if (times[[i]][l] <= setpoint1) c2[l] <- 1
        else if (times[[i]][l] <= setpoint2) c2[l] <- 1 + 1*sin(1*pi*(times[[i]][l] - setpoint1)/(setpoint2 - setpoint1))
        else if (times[[i]][l] <= setpoint3) c2[l] <- 1
        else if (times[[i]][l] <= setpoint4) c2[l] <- 1
        else c2[l] <- c2[l] <- 1
      }
      epsilon <- mvrnorm(1, mu = mu, Sigma = sig)
      oridata_list[[i]] <- c2 + epsilon
    }

    if (i %in% class3) {
      c3 <- 0*times[[i]]
      for (l in 1:n0) {
        if (times[[i]][l] <= setpoint1) c3[l] <- 1
        else if (times[[i]][l] <= setpoint2) c3[l] <- 1 - 1*sin(1*pi*(times[[i]][l] - setpoint1)/(setpoint2 - setpoint1))
        else if (times[[i]][l] <= setpoint3) c3[l] <- 1
        else if (times[[i]][l] <= setpoint4) c3[l] <- 1 + 1*sin(1*pi*(times[[i]][l] - setpoint3)/(setpoint4 - setpoint3))
        else c3[l] <- 1
      }
      epsilon <- mvrnorm(1, mu = mu, Sigma = sig)
      oridata_list[[i]] <- c3 + epsilon
    }


    if (i %in% class4) {
      c4 <- 0*times[[i]]
      for (l in 1:n0) {
        if (times[[i]][l] <= setpoint1) c4[l] <- 1
        else if (times[[i]][l] <= setpoint2) c4[l] <- 1 - 1*sin(1*pi*(times[[i]][l] - setpoint1)/(setpoint2 - setpoint1))
        else if (times[[i]][l] <= setpoint3) c4[l] <- 1 + 1*sin(1*pi*(times[[i]][l] - setpoint2)/(setpoint3 - setpoint2))
        else if (times[[i]][l] <= setpoint4) c4[l] <- 1 - 1*sin(1*pi*(times[[i]][l] - setpoint3)/(setpoint4 - setpoint3))
        else c4[l] <- 1
      }

      epsilon <- mvrnorm(1, mu = mu, Sigma = sig)
      oridata_list[[i]] <- c4 + epsilon
    }
  }

  return(list(location = times, data_lst = oridata_list, clus.true = c(class1, class2, class3, class4)))
}
