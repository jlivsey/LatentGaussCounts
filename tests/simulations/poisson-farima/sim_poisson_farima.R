n.seq   = c(100, 200, 400)
d.seq = c(-.75, -.25, .25, .75)
lam.seq = c(2, 5, 10)
Nsim = 1000
total.iter = Nsim * length(n.seq) * length(d.seq) * length(lam.seq)
save.cols = c("estim.method", "n",
              "lam.true", "lam.est", "lam.se",
              "d.true", "d.est", "d.se")
simResults_poisson_farima = data.frame(matrix(nrow = total.iter,
                                              ncol = length(save.cols)))
colnames(simResults.mixedpoisAR1) = save.cols

idx = 1
pb <- txtProgressBar(min=2,max=total.iter,style=3)
for(n.idx in 1:length(n.seq)){
  for(phi.idx in 1:length(d.seq)){
    for(lam.idx in 1:length(lam.seq)){
      for(dumb.variable in 1:Nsim){
        n   = n.seq[n.idx]
        phi = d.seq[phi.idx]
        lam = lam.seq[lam.idx]
        x   = sim_pois_ar(n = n, phi = phi, lam = lam)
        optim.output = LGC(x,
                        count.family = "Poisson",
                        gauss.series = "FARIMA",
                        estim.method = "gaussianLik")
        simResults_poisson_farima$estim.method[idx] = "gaussianLik"
        simResults_poisson_farima$n[idx] = n
        simResults_poisson_farima$lam.true[idx] = lam
        simResults_poisson_farima$lam.est[idx]  = optim.output$par[1]
        simResults.mixedpoisAR1$lam.se[idx]   = optim.output$stder[1]
        simResults_poisson_farima$d.true[idx] = phi
        simResults_poisson_farima$d.est[idx]  = optim.output$par[2]
        simResults_poisson_farima$d.se[idx]   = optim.output$stder[2]
        idx = idx + 1
        setTxtProgressBar(pb,idx)
      }}}}
close(pb)

save(simResults.poisAR1, file = "simResults_poisAR1.Rdata")




# ---- Windows parallel -------------------------------------------------------
library(parallel)
# load data as list
l2 <- list()
for(i in 1:100)     l2[[i]] = sim_pois_ar(n = 200, phi = 0.7,
                                          lam = 2)
for(i in 101:200)  l2[[i]] = sim_pois_ar(n = 200, phi = 0.2,
                                         lam = 2)
for(i in 201:300) l2[[i]] = sim_pois_ar(n = 200, phi = 0.7,
                                        lam = 10)
for(i in 301:400) l2[[i]] = sim_pois_ar(n = 200, phi = -0.7,
                                        lam = 2)
# set up cluster
cl <- makeCluster(detectCores())
# load our package for each node
clusterEvalQ(cl, library(latentGaussCounts))
# export data to each node
clusterExport(cl, varlist = "l2")
# Run LGC()
system.time({
  out <- parLapply(cl, l2[1:5], function(x){
    LGC(x,
        count.family = "Poisson",
        gauss.series = "AR", p=1,
        estim.method = "particlesSIS")
  })
})
# First 5
system.time({
  out2 <- parLapply(cl, l2[6:25], function(x){
    LGC(x,
        count.family = "Poisson",
        gauss.series = "AR", p=1,
        estim.method = "particlesSIS")
  })
})
# Next 20
system.time({
  out3 <- parLapply(cl, l2[26:100], function(x){
    LGC(x,
        count.family = "Poisson",
        gauss.series = "AR", p=1,
        estim.method = "particlesSIS")
  })
})
system.time({
  out4 <- parLapply(cl, l2[101:200], function(x){
    LGC(x,
        count.family = "Poisson",
        gauss.series = "AR", p=1,
        estim.method = "particlesSIS")
  })
})
system.time({
  out5 <- parLapply(cl, l2[201:300], function(x){
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

