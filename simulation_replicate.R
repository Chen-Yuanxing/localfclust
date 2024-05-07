
###############################################################################################################
##############    This file is the main file for reproducing simulation results            ####################
##############    of "Local clustering for functional data" by Chen et al. (2024).         ####################
##############    Here we take replicating the results of Example 1 as an example,         ####################
##############    and similar procedures for other examples can be obtained accordingly.      ####################
###############################################################################################################


# Example 1
library(localfclust)

data.example1 <- SimulateExample1(seed = 1)  # generate simulated data for Example 1


# determine the optimal lambda1
lambda1.candidate <- c(seq(1e-7, 1e-3, length.out = 100))
lambda1.opt <- CalculateGCV(times = data.example1$location, oridata = data.example1$data_lst, 
                            lambda1_seq = lambda1.candidate)  



# determine the optimal lambda2 and lambda3
lambda2.candidate <- seq(1, 2, 1)
lambda3.candidate <- seq(0.05, 0.5, 0.05)
bic.result <- CalculateBIC(times = data.example1$location, oridata = data.example1$data_lst, lambda1 = lambda1.opt, 
                           lambda2_seq = lambda2.candidate, lambda3_seq = lambda3.candidate)  


nlambda2 <- length(lambda2.candidate)
nlambda3 <- length(lambda3.candidate)
bic.mat <- matrix(0, nrow = nlambda2, ncol = nlambda3)
for (l1 in 1:nlambda2) {
  for (l2 in 1:nlambda3) {
    bic.mat[l1, l2] <- bic.result[l1, l2][[1]]$bic_val
  }
}

id <- which(bic.mat == min(bic.mat), arr.ind = T)[1, ]
lambda2.opt <- lambda2.candidate[id[1]]
lambda3.opt <- lambda3.candidate[id[2]]


# model fitting process
model.lfc <- LFC.fit(times = data.example1$location, oridata = data.example1$data_lst, lambda1 = lambda1.opt, 
                     lambda2 = lambda2.opt, lambda3 = lambda3.opt)




plot.individual.localfclust(model.lfc = model.lfc)  # show the figure of individual smoothed curves.
plot.center.localfclust(model.lfc) # show the figure of mean curves for obtained cluster.





# result evaluation process
library(flexclust)
library(funtimes)
library(fda)
library(igraph)

EvaluateResultsForExample1 <- function(setpoint1 = 0.6, n = 100, n0 = 100, model.lfc, data.example1){
  
  N <- n0*n
  
  basisobj <- model.lfc$basisobj
  p <- nrow(model.lfc$Beta)
  
  Z_lst <- vector(mode = "list", length = n)
  for (i in 1:n) {
    Z_lst[[i]] <- eval.basis(data.example1$location[[i]], basisobj)
  }
  
  class1 <- 1:(n/2)
  class2 <- (n/2 + 1):n
  
  setpoint2 <- 1
  
  
  c1 <- rep(0, n0)
  for (l in 1:n0) {
    if (Times[l] >= setpoint1 & Times[l] <= setpoint2) c1[l] <- 1 + 1*sin(2*pi*(Times[l] - setpoint1)/(setpoint2 - setpoint1))
    else c1[l] <- 1
  }
  
  c2 <- rep(0, n0)
  for (l in 1:n0) {
    if (Times[l] >= setpoint1 & Times[l] <= setpoint2) c2[l] <- 1 - 1*sin(2*pi*(Times[l] - setpoint1)/(setpoint2 - setpoint1))
    else c2[l] <- 1 
  }
  
  Beta <- model.lfc$Beta
  Alpha <- model.lfc$centers
  
  Khat <- ncol(Alpha)
  
  group <- rep(0, n)
  for (i in 1:n) {
    dis_vec <- rep(0, Khat)
    for (k in 1:Khat) {
      dis_vec[k] <- sum((Beta[,i] - Alpha[,k])^2)
    }
    group[i] <- which.min(dis_vec)
  }
  

  Rand_index <- comPart(group, c(rep(1, n/2), rep(2, n/2)))[2]
  Purity_index <- purity(group, c(rep(1, n/2), rep(2, n/2)))[[1]]
  
  
  mse <- rep(0, n)
  for (i in 1:n) {
    if (i %in% class1)  mse[i] <- sum((Z_lst[[i]] %*% Beta[,i] - c1)^2)
    else  mse[i] <- sum((Z_lst[[i]] %*% Beta[,i] - c2)^2)
  }
  RMSE <- sqrt(sum(mse)/N)
  
  
  
  aa <- seq(0, 1, length.out = n0)
  LL <- eval.basis(aa, basisobj)
  
  point_mat <- matrix(0, nrow = n0, ncol = n)
  for (i in 1:n) {
    point_mat[,i] <- as.numeric(LL %*% Beta[,i])
  }
  
  #interval_Khat <- rep(0, 100)
  interval_K <- rep(0, n0)
  interval_RI <- rep(0, n0)
  interval_JI <- rep(0, n0)
  interval_PF <- rep(0, n0)
  membership_lst <- vector(mode = "list", length = n0)
  Mhat <- 1
  
  #average_ari <- rep(0, 100)
  for (j in 1:n0) {
    dist_mat <- as.matrix(dist(point_mat[j,]))
    adj_vec <- dist_mat[t(upper.tri(dist_mat))]
    adj_mat <- matrix(adj_vec, nrow = 1)
    Ad_final <- CreateAdjacency(adj_mat, n);
    G_final <- graph.adjacency(Ad_final, mode = 'upper')
    cls_final <- components(G_final)
    interval_K[j] <- cls_final$no
    membership_lst[[j]] <- cls_final$membership
    
    if (j > (n0*0.05) & j < (n0*0.95)){
      if (purity(membership_lst[[j]], membership_lst[[j-1]])[[1]] != 1) Mhat <- Mhat + 1
    }
    
    for (k in 1:cls_final$no) {
      point_mat[j, cls_final$membership == k] <- k
    }
    if (aa[j] <= setpoint1) {
      interval_RI[j] <- comPart(point_mat[j,], rep(1,n))[2]
      interval_PF[j] <- purity(point_mat[j,], rep(1,n))[[1]]
    } 
    
    else {
      interval_RI[j] <- comPart(point_mat[j,], c(rep(1,n/2), rep(2, n/2)))[2]
      interval_PF[j] <- purity(point_mat[j,], c(rep(1,n/2), rep(2, n/2)))[[1]]
    }  
    
  }
  
  K1hat <- mean(interval_K[1+(n0*0.05):(n0*setpoint1)])
  RI_local1 <- mean(interval_RI[1+(n0*0.05):(n0*setpoint1)])
  PF_local1 <- mean(interval_PF[1+(n0*0.05):(n0*setpoint1)])
  
  K2hat <- mean(interval_K[(n0*setpoint1+1):(n0*0.95)])
  RI_local2 <- mean(interval_RI[(n0*setpoint1+1):(n0*0.95)])
  PF_local2 <- mean(interval_PF[(n0*setpoint1+1):(n0*0.95)])
  
  RI_local_average <- mean(interval_RI[1+(n0*0.05):(n0*0.95)])
  PF_local_average <- mean(interval_PF[1+(n0*0.05):(n0*0.95)])
  

  
  
  index_list = list(RMSE_LFC = RMSE, Khat_LFC = Khat, Rand_LFC = Rand_index, Purity_LFC = Purity_index, 
                    RI_local1_LFC = RI_local1, PF_local1_LFC = PF_local1, K1hat_LFC = K1hat,
                    RI_local2_LFC = RI_local2, PF_local2_LFC = PF_local2, K2hat_LFC = K2hat, Mhat_LFC = Mhat,
                    RI_local_average_LFC = RI_local_average, PF_local_average_LFC = PF_local_average)
  
  return(index_list)
}

evaluation.result <- EvaluateResultsForExample1(model.lfc = model.lfc, data.example1 = data.example1)

