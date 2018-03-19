#set.seed(Sys.time())
x=sim_pois_ar(n = 200, phi = 0.7, lam = 2)
plot.ts(x)
acf(x,lag.max = 20)

LGC(x, count.family = "Poisson",
       gauss.series = "AR", p=1,
       estim.method = "gaussianLik")

LGC(x, count.family = "Poisson",
       gauss.series = "AR", p=1,
       estim.method = "particlesSIS")

