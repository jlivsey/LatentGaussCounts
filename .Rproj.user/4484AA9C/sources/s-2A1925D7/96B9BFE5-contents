
sim_mixedpois_ar1 = function(n, phi, p, lam1, lam2){
  z = arima.sim(model = list(ar=phi), n = n); z = z/sd(z) # standardized
  x = qmixpois(y = pnorm(z), p = p, lam1 = lam1, lam2 = lam2)
  return(x)
}


x = sim_mixedpois_ar1(200, .7, .5, 1, 15)
LGC(x,
    count.family = "mixed-Poisson", n.mix=2,
    gauss.series = "AR", p=1,
    estim.method = "gaussianLik",
    print.progress = TRUE)

LGC(x,
    count.family = "mixed-Poisson", n.mix=2,
    gauss.series = "AR", p=1,
    estim.method = "particlesSIS",
    print.progress = TRUE)
