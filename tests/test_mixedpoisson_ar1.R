
x = sim_mixedpois_ar1(n = 400, phi = .7, p = .5, lam1 = 2, lam2 = 10)

LGC(x,
    count.family = "mixed-Poisson", n.mix=2,
    gauss.series = "AR", p=1,
    estim.method = "gaussianLik",
    print.progress = TRUE, print.initial.estimates = TRUE)

LGC(x,
    count.family = "mixed-Poisson", n.mix=2,
    gauss.series = "AR", p=1,
    estim.method = "particlesSIS",
    print.progress = TRUE)
