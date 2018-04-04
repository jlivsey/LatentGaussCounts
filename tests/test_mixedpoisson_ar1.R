
x = sim_mixedpois_ar1(400, .7, .5, 1, 15)

LGC(x,
    count.family = "mixed-Poisson", n.mix=2,
    gauss.series = "AR", p=1,
    estim.method = "gaussianLik",
    print.progress = FALSE, print.initial.estimates = TRUE)

LGC(x,
    count.family = "mixed-Poisson", n.mix=2,
    gauss.series = "AR", p=1,
    estim.method = "particlesSIS",
    print.progress = TRUE)


