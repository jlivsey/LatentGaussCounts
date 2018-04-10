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

# ---- parallel test ----------------------------------------------------------

library(parallel)
# load data as list
l2 <- list()
l2[[1]] = x
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
                                          count.family = "Poisson",
                                          gauss.series = "AR", p=1,
                                          estim.method = "particlesSIS")
                                    })
})
# Stop the cluster
stopCluster(cl)
