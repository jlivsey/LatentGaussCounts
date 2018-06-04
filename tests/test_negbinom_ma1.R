
n = 400
theta = .25
# theta = c(0.1,0.7)
r = 3
p = 1/2
x = sim_negbinom_ma1(n, theta, r, p)
# plot.ts(x)
# acf(x,lag.max = 20)
# pacf(x, lag.max=20)

LGC(x,
    count.family = "negbinom",
    gauss.series = "MA", q=1,
    estim.method = "gaussianLik",
    print.initial.estimates = TRUE,
    print.progress = TRUE)

# ---- quick sim ----

n = 400
theta = .25
r = 3
p = 1/2

Nsim = 5
params = matrix(ncol=3, nrow=Nsim)
ses = matrix(ncol=3, nrow=Nsim)
for(i in 1:Nsim)
{
x <<- sim_negbinom_ma1(n, theta, r, p)
optim.out = LGC(x,
    count.family = "negbinom",
    gauss.series = "MA", q=1,
    estim.method = "gaussianLik",
    print.initial.estimates = TRUE,
    print.progress = FALSE)
params[i,] = optim.out$par
ses[i,] = optim.out$stder
}
