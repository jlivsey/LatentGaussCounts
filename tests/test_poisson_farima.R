x = pois.farima(Tt = 200, d = .4, lam = 2)
plot.ts(x)
acf(x,lag.max = 50)


LGC(x,
    count.family = "Poisson",
    gauss.series = "FARIMA",
    estim.method = "gaussianLik")

LGC(x,
    count.family = "Poisson",
    gauss.series = "FARIMA",
    estim.method = "particlesSIS")
