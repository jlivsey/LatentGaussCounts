# Brockwell and Davis UG book, pp. 90-93
xi1 = 2
xi2 = 5

xi1 = 2*(1+1i*sqrt(3))/3
xi2 = 2*(1-1i*sqrt(3))/3

phi1 = Re(1/xi1 + 1/xi2)
phi2 = -Re((1/xi1) * (1/xi2))

c(phi1,phi2)

#set.seed(Sys.time())
x=sim_pois_ar1(n = 200, phi = c(phi1,phi2), lam = 2)
plot.ts(x)
acf(x,  lag.max = 20)
pacf(x, lag.max = 20)

LGC(x, count.family = "Poisson",
    gauss.series = "AR", p=2,
    estim.method = "gaussianLik", print.progress = TRUE)

LGC(x, count.family = "Poisson",
    gauss.series = "AR", p=2,
    estim.method = "particlesSIS")

