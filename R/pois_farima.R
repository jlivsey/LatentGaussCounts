#' Simulate from Poisson-FARIMA(0,d,0) model
#'
#' @param Tt length of series to be generated
#' @param d long-memory parameter
#' @param lam poisson parameter
#'
#' @return discrete valued time series with Poisson marginals and FARIMA(0,d,0)
#'     covariance structure
#' @export
#'
#' @examples pois.farima(200, 1/4, 3)


pois.farima <- function(Tt,d,lam){
  zt <- farima.sim(Tt,d)*gamma(1-d)/sqrt(gamma(1-2*d))
  xt <- qpois(pnorm(zt,0,1),lam)
  return(xt)
}
