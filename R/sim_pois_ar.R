#' Simulate from Poisson-AR1 model
#'
#' @param n length of series
#' @param phi AR(1) parameter
#' @param lam Poisson parameter
#'
#' @return simulated count series of length n
#' @export
#'

sim_pois_ar = function(n, phi, lam){
  z = arima.sim(model = list(ar=phi), n = n); z = z/sd(z) # standardized
  x = qpois(pnorm(z), lam)
  return(x)
}
