n.seq   = c(100, 200, 400)
phi.seq = c(-.75, -.25, .25, .75)
lam.seq = c(2, 5, 10)
Nsim = 1000
total.iter = Nsim * length(n.seq) * length(phi.seq) * length(lam.seq)
save.cols = c("estim.method", "n",
              "lam.true", "lam.est", "lam.se",
              "phi.true", "phi.est", "phi.se")
simResults.poisAR1 = matrix(nrow = total.iter, ncol = length(save.cols))
colnames(simResults.poisAR1) = save.cols

idx = 1
pb <- txtProgressBar(min=2,max=total.iter,style=3)
for(n.idx in 1:length(n.seq)){
for(phi.idx in 1:length(phi.seq)){
for(lam.idx in 1:length(lam.seq)){
for(dumb.variable in 1:Nsim){
  n   = n.seq[n.idx]
  phi = phi.seq[phi.idx]
  lam = lam.seq[lam.idx]
  x   = sim_pois_ar(n = n, phi = phi, lam = lam)
  param.est = LGC(x,
                  count.family = "Poisson",
                  gauss.series = "AR", p=1,
                  estim.method = "gaussianLik")$par
  simResults.poisAR1[idx, "estim.method"] = "gaussianLik"
  simResults.poisAR1[idx, "n"] = n
  simResults.poisAR1[idx, "lam.true"] = lam
  simResults.poisAR1[idx, "lam.est"]  = param.est[1]
  simResults.poisAR1[idx, "phi.true"] = phi
  simResults.poisAR1[idx, "phi.est"]  = param.est[2]
  idx = idx + 1
  setTxtProgressBar(pb,idx)
}}}}
close(pb)

save(simResults.poisAR1, file = "simResults_poisAR1.Rdata")




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

