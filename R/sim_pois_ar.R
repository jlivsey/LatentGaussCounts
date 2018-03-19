

sim_pois_ar = function(n, phi, lam){
  z = arima.sim(model = list(ar=phi), n = n); z = z/sd(z) # standardized
  x = qpois(pnorm(z), lam)
  return(x)
}
