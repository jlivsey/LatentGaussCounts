# library(orthopolynom) # for hermite polynomials
library(mvtnorm) # for multivariate Gaussian density

sim_pois_ar1 = function(n, phi, lam){
  z = arima.sim(model = list(ar=phi), n = n); z = z/sd(z) # standardized
  x = qpois(pnorm(z), lam)
  return(x)
}

x=sim_pois_ar1(200,0.7,2)
LGC(x,p=1,estim.method = "gaussianLik")
LGC(x,p=1,estim.method = "particlesSIS")


v = replicate(n = 20, expr = {
  x = sim_pois_ar1(200, .7, 2)
  LGC(x, p=1)$par
})

