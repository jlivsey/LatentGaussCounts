#' Simulate from Poisson-MA model
#'
#' @param n length of series
#' @param theta MA parameter
#' @param lam Poisson parameter
#'
#' @return simulated count series of length n
#' @export
#'

sim_pois_ma = function(n, theta, lam){
  z = arima.sim(model = list(ma=theta), n = n); z = z/sd(z) # standardized
  x = qpois(pnorm(z), lam)
  return(x)
}
