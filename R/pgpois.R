#' cdf of generalized poisson
#'
#' @param x series
#' @param lambda 1 parameter
#' @param eta 2 parameter
#'
#' @return cdf
#' @export
#'

pgpois= function(x,lambda, eta) { # using Yisu's code for the cdf function
  cdf.vec <- rep(-99,length(x))
  for (i in 1:length(x)){
    cdf.vec[i] <- sum(dgpois(0:x[i],lambda,eta))
  }
  return(cdf.vec)
}


