
n = 400
theta = .9
r = 3
p = 1/2
x = sim_negbinom_ma1(n, theta, r, p)
plot.ts(x)
acf(x,lag.max = 20)
pacf(x, lag.max=20)

LGC(x,
    count.family = "negbinom",
    gauss.series = "MA", q=1,
    estim.method = "gaussianLik",
    print.initial.estimates = TRUE,
    print.progress = TRUE)
