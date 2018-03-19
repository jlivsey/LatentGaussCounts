# library(orthopolynom) # for hermite polynomials
# library(mvtnorm) # for multivariate Gaussian density

sim_pois_ar1 = function(n, phi, lam){
  z = arima.sim(model = list(ar=phi), n = n); z = z/sd(z) # standardized
  x = qpois(pnorm(z), lam)
  return(x)
}

#set.seed(Sys.time())
x=sim_pois_ar1(200,0.7,2)
plot.ts(x)
acf(x,lag.max = 20)

LGC(x, count.family = "Poisson",
       gauss.series = "AR", p=1,
       estim.method = "gaussianLik")

LGC(x, count.family = "Poisson",
       gauss.series = "AR", p=1,
       estim.method = "particlesSIS")


# Brockwell and Davis UG book, pp. 90-93
xi1 = 2
xi2 = 5

xi1 = 2*(1+1i*sqrt(3))/3
xi2 = 2*(1-1i*sqrt(3))/3

phi1 = Re(1/xi1 + 1/xi2)
phi2 = -Re((1/xi1) * (1/xi2))

c(phi1,phi2)

#set.seed(Sys.time())
x=sim_pois_ar1(200,c(phi1,phi2),2)
plot.ts(x)
acf(x,lag.max = 20)

LGC(x, count.family = "Poisson",
    gauss.series = "AR", p=2,
    estim.method = "particlesSIS")

