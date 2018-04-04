#' Simulate from mixed-Poisson-AR1 model
#'
#' @param n number of observations
#' @param phi AR parameter
#' @param p mixing probability
#' @param lam1 mean of first Poisson component
#' @param lam2 mean of second Poisson component
#'
#' @return simulated mixed-Poisson-AR1 model of length n
#' @export
#'


sim_mixedpois_ar1 = function(n, phi, p, lam1, lam2){
  z = arima.sim(model = list(ar=phi), n = n); z = z/sd(z) # standardized
  x = qmixpois(y = pnorm(z), p = p, lam1 = lam1, lam2 = lam2)
  return(x)
}
