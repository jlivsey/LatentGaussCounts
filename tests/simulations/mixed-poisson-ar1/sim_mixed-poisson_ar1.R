n.seq    = c(200, 400)
phi.seq  = c(.25, .75)
p.seq    = 1/4
lam1.seq = 2
lam2.seq = c(3, 5, 10)
Nsim = 1000
total.iter = Nsim * length(n.seq) * length(phi.seq) *
             length(lam1.seq) * length(lam2.seq) * length(p.seq)
save.cols = c("estim.method", "n",
              "lam1.true", "lam1.est", "lam1.se",
              "lam2.true", "lam2.est", "lam2.se",
              "p.true", "p.est", "p.se",
              "phi.true", "phi.est", "phi.se")
simResults.mixedpoisAR1 = data.frame(matrix(nrow = total.iter,
                                            ncol = length(save.cols)))
colnames(simResults.mixedpoisAR1) = save.cols

idx = 1
pb <- txtProgressBar(min=2,max=total.iter,style=3)
for(n.idx in 1:length(n.seq)){
for(phi.idx in 1:length(phi.seq)){
for(lam1.idx in 1:length(lam1.seq)){
for(lam2.idx in 1:length(lam2.seq)){
for(p.idx in 1:length(p.seq)){
for(dumb.variable in 1:Nsim){
  n   = n.seq[n.idx]
  phi = phi.seq[phi.idx]
  lam1 = lam1.seq[lam1.idx]
  lam2 = lam2.seq[lam2.idx]
  p = p.seq[p.idx]
  x = sim_mixedpois_ar1(n = n, phi = phi, p = p, lam1 = lam1, lam2 = lam2)
  optim.output = LGC(x,
                     count.family = "mixed-Poisson", n.mix=2,
                     gauss.series = "AR", p=1,
                     estim.method = "gaussianLik")
  simResults.mixedpoisAR1$estim.method[idx] = "gaussianLik"
  simResults.mixedpoisAR1$n[idx]         = n
  simResults.mixedpoisAR1$p.true[idx]    = p
  simResults.mixedpoisAR1$p.est[idx]     = optim.output$par[1]
  simResults.mixedpoisAR1$p.se[idx]      = optim.output$stder[1]
  simResults.mixedpoisAR1$lam1.true[idx] = lam1
  simResults.mixedpoisAR1$lam1.est[idx]  = optim.output$par[2]
  simResults.mixedpoisAR1$lam1.se[idx]   = optim.output$stder[2]
  simResults.mixedpoisAR1$lam2.true[idx] = lam2
  simResults.mixedpoisAR1$lam2.est[idx]  = optim.output$par[3]
  simResults.mixedpoisAR1$lam2.se[idx]   = optim.output$stder[3]
  simResults.mixedpoisAR1$phi.true[idx]  = phi
  simResults.mixedpoisAR1$phi.est[idx]   = optim.output$par[4]
  simResults.mixedpoisAR1$phi.se[idx]    = optim.output$stder[4]
  idx = idx + 1
  setTxtProgressBar(pb,idx)
}}}}}}
close(pb)


# ---- Windows parallel -------------------------------------------------------
library(parallel)
# load data as list
l2 <- list()
for(i in 1:500) l2[[i]] = sim_mixedpois_ar1(n = 200, phi = .7, p = .25, lam1 = 2, lam2 = 3)
for(i in 501:1000) l2[[i]] = sim_mixedpois_ar1(n = 200, phi = .7, p = .25, lam1 = 2, lam2 = 10)
for(i in 1001:1500) l2[[i]] = sim_mixedpois_ar1(n = 200, phi = .2, p = .25, lam1 = 2, lam2 = 10)
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

