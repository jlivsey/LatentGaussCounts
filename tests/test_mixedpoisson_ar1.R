
x = sim_mixedpois_ar1(n = 200, phi = .7, p = .5, lam1 = 2, lam2 = 10)

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

# ---- parallel test ----------------------------------------------------------

library(parallel)
# load data as list
l2 <- list()
l2[[1]] <- x
# set up cluster
cl <- makeCluster(detectCores())
# load our package for each node
clusterEvalQ(cl, library(latentGaussCounts))
# export data to each node
clusterExport(cl, varlist = "l2")
# Run LGC()
system.time({
out <- parLapply(cl, l2, function(x){
                                      LGC(x,
                                          count.family = "mixed-Poisson", n.mix=2,
                                          gauss.series = "AR", p=1,
                                          estim.method = "particlesSIS")
                                    })
})
# Stop the cluster
stopCluster(cl)
