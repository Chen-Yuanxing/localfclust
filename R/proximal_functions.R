#' @title Proximal mapping functions
#' @description The proximal (average) mapping of mcp function.
#'
#'
#' @param x The input vector.
#' @param l An integer taking values from \{1,...,L\} and L is the number of subintervals.
#' @param r An integer specifying the degree of B-splines.
#' @param varsigma A positive double scaling parameter.
#' @param tau A positive double parameter that controls the concavity of mcp function. The default value is 3.
#' @param lambda The tuning parameter for mcp function.
#'
#' @return A vector with the same dimensions as \code{x}.
#' @export ProximalMCP
#' @export ProximalAverageMCP
#'


#' @name proximal
NULL

#' @rdname proximal
ProximalMCP <- function(x, varsigma, tau = 3, lambda){

  p <- length(x)
  x_norm <- sqrt(sum(x^2))
  if (x_norm <= lambda*varsigma) z <- rep(0, p)
  else if (x_norm <= lambda*tau) z <- max(0, x_norm - lambda*varsigma)/(1 - varsigma/tau) * x/x_norm
  else z <- x
}

#' @rdname proximal
ProximalAverageMCP <- function(x, l, r, varsigma, tau = 3, lambda){

  p <- length(x)
  xl_norm <- sqrt(sum(x[l:(l+r)]^2))
  z <- rep(0, p)

  for (j in 1:p) {
    if ((j %in% l:(l+r)) & (xl_norm <= lambda*varsigma)) z[j] <- 0
    if ((j %in% l:(l+r)) & (xl_norm <= lambda*tau)) z[j] <- tau/(tau-varsigma) * max((1-varsigma*lambda/xl_norm), 0) * x[j]
    else z[j] <- x[j]
  }

  return(z)
}



