#set.seed(Sys.time())
x=sim_pois_ma(n = 200, theta = 0.5, lam = 2)
#plot.ts(x)
#acf(x,lag.max = 20)

LGC(x,
    count.family = "Poisson",
    gauss.series = "MA", q=1,
    estim.method = "gaussianLik",
    print.initial.estimates = TRUE,
    print.progress = TRUE)

# LGC(x, count.family = "Poisson",
#     gauss.series = "AR", p=1,
#     estim.method = "particlesSIS")
