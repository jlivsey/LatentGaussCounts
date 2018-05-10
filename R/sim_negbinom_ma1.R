#' Simulate from NegativeBinomial-MA model
#'
#' @param n length of series
#' @param theta MA(1) parameter
#' @param r number of failures before x successes
#' @param p success probability
#'
#' @return simulated count series of length n
#' @export
#'

sim_negbinom_ma1 = function(n, theta, r, p){
  z = arima.sim(model = list(ma=theta), n = n); z = z/sd(z) # standardized
  x = qnbinom(p = pnorm(z), size = r, prob = p)
  return(x)
}
