#set.seed(Sys.time())
x=sim_pois_ar(n = 400, phi = 0.7, lam = 2)
plot.ts(x)
acf(x,lag.max = 20)

LGC(x,
    count.family = "Poisson",
    gauss.series = "AR", p=1,
    estim.method = "gaussianLik",
    print.initial.estimates = TRUE,
    print.progress = TRUE)

LGC(x, count.family = "Poisson",
       gauss.series = "AR", p=1,
       estim.method = "particlesSIS")

