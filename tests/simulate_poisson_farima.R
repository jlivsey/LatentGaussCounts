
xt = pois.farima(200,0.4,2)
plot.ts(xt)
acf(xt)


LGC(x, count.family = "Poisson",
       gauss.series = "FARIMA",
       estim.method = "particlesSIS")


