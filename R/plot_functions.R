#' @title Plot the results of  the LFC method
#' @description Provides plots of the classified curves and estimated cluster mean functions, respectively.
#'
#' @param model.lfc The output of \code{\link{LFC.fit}}.

#' @return The ggplot2 object.

#'
#' @import fda
#' @import ggplot2
#' @export plot.individual.localfclust
#' @export plot.center.localfclust
#'
#' @name Plot.LFC
NULL

#' @rdname Plot.LFC
plot.individual.localfclust <- function(model.lfc) {
  Beta_ini <- model.lfc$Beta_ini
  basis.fd <- model.lfc$basisobj

  n <- ncol(Beta_ini)
  n0 <- 1000
  Times <- seq(0, 1, len = n0)

  Z1 <- eval.basis(Times, basis.fd)

  y.individual.hat.mat <- matrix(0, nrow = n0, ncol = n)
  for (i in 1:n) {
    y.individual.hat.mat[, i] <- as.numeric(Z1 %*% Beta_ini[, i])
  }

  fit.individual.data <- data.frame(
    Times.all = rep(Times, n),
    value.all = c(as.vector(y.individual.hat.mat)),
    Individuals = as.character(c(rep(1:n, each = n0)))

  )
  fit.individual.data$Individuals <- factor(fit.individual.data$Individuals, levels = as.character(1:n))

  plot.individual <- ggplot(fit.individual.data, aes(x = Times.all, y = value.all, color = Individuals, linetype = Individuals)) +
    geom_line()  +
    scale_color_manual(values = model.lfc$cls_mem) +
    scale_linetype_manual(values = model.lfc$cls_mem) +
    xlab("Time") +
    ylab("") +
    theme(axis.title.y = element_text(size = 15), axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10), axis.title.x = element_text(size = 15))

  return(plot.individual)
}

#' @rdname Plot.LFC
plot.center.localfclust <- function(model.lfc) {

  Alpha <- model.lfc$centers
  basis.fd <- model.lfc$basisobj


  n0 <- 1000
  Times <- seq(0, 1, len = n0)

  Z1 <- eval.basis(Times, basis.fd)
  K.hat <- model.lfc$cls_num

  y.center.hat.mat <- matrix(0, nrow = n0, ncol = K.hat)
  for (k in 1:K.hat) {
    y.center.hat.mat[, k] <- as.numeric(Z1 %*% Alpha[, k])
  }

  fit.center.data <- data.frame(
    Times.all = rep(Times, K.hat),
    value.all = c(as.vector(y.center.hat.mat)),
    Clusters = as.character(c(rep(1:K.hat, each = n0)))
  )
  fit.center.data$Clusters <- factor(fit.center.data$Clusters, levels = as.character(1:K.hat), labels = paste0("Cluster ", 1:K.hat))

  plot.center <- ggplot(fit.center.data, aes(x = Times.all, y = value.all, color = Clusters, linetype = Clusters)) +
    geom_line()  +
    xlab("Time") +
    ylab("") +
    theme(axis.title.y = element_text(size = 15), axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10), axis.title.x = element_text(size = 15))

  return(plot.center)
}

