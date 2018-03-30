#set.seed(Sys.time())
x=sim_pois_ar(n = 200, phi = 0.7, lam = 2)
plot.ts(x)
acf(x,lag.max = 20)

LGC(x,
    count.family = "Poisson",
    gauss.series = "AR", p=1,
    estim.method = "gaussianLik")

LGC(x, count.family = "Poisson",
       gauss.series = "AR", p=1,
       estim.method = "particlesSIS")


param.est = matrix(nrow=50, ncol=2)
colnames(param.est) <- c("lambda", "phi")
for(i in 1:50){
  if(i%%10==0) print(i)
  x=sim_pois_ar(n = 200, phi = 0.2, lam = 2)
  param.est[i, ] = LGC(x,
                       count.family = "Poisson",
                       gauss.series = "AR", p=1,
                       estim.method = "gaussianLik")$par
}

{
op = par(mfrow=c(1,2), mar=c(1, 3, 3, 1))
boxplot(param.est[,1], main="lambda (true value 2)"); abline(h=2, lty="dotted", col=2)
boxplot(param.est[,2], main="phi (true value 0.2)"); abline(h=0.2, lty="dotted", col=2)
par(op)
}


# ---- Windows parallel -------------------------------------------------------
library(parallel)
# load data as list
l2 <- list()
for(i in 1:10) l2[[i]] = sim_pois_ar(n = 200, phi = 0.7, lam = 2)
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


# ---- Linux or OS-X parallel -------------------------------------------------------
library(parallel)
# load data as list
l2 <- list()
for(i in 1:10) l2[[i]] = sim_pois_ar(n = 200, phi = 0.7, lam = 2)
# Run LGC() function
mclapply(l2, function(x){
                          LGC(x,
                              count.family = "Poisson",
                              gauss.series = "AR", p=1,
                              estim.method = "particlesSIS")
                        }, mc.cores = detectCores())

