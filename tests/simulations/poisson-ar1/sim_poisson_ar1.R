n.seq   = c(100, 200, 400)
phi.seq = c(-.75, -.25, .25, .75)
lam.seq = c(2, 5, 10)
Nsim = 200
total.iter = Nsim * length(n.seq) * length(phi.seq) * length(lam.seq)
save.cols = c("estim.method", "n",
              "lam.true", "lam.est", "lam.se",
              "phi.true", "phi.est", "phi.se")
simResults.poisAR1 = matrix(nrow = total.iter, ncol = length(save.cols))
colnames(simResults.poisAR1) = save.cols


# --- Add Gaussian Likelihood results to simResults matrix ----
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

# save(simResults.poisAR1, file = "simResults_poisAR1.Rdata")



# ---- Add Particle filtering results to SimResults matrix ----
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


# ---- Add implied Yule-Walker estim results to SimResults matrix ----
library(latentGaussCounts)
library(orthopolynom)
polys <- hermite.he.polynomials(100) # In global env for g_coefs() function - BAD PROGRAMING
source("~/Dropbox/jim/MVcopula/code/functions.R")
n.coefs = 20

# Current size of previous simulation results
idx <- dim(simResults_poisson_ar1)[1] + 1

# Set current results DF equal to old restuls
simResults.poisAR1 <- simResults_poisson_ar1
# Coerce 'estim.method' from factor to character
simResults.poisAR1$estim.method <- as.character(simResults.poisAR1$estim.method)

pb <- txtProgressBar(min = idx,
                     max = total.iter + idx,
                     style = 3)
for(n.idx in 1:length(n.seq)){
  for(phi.idx in 1:length(phi.seq)){
    for(lam.idx in 1:length(lam.seq)){
      for(dumb.variable in 1:Nsim){
        n   = n.seq[n.idx]
        phi = phi.seq[phi.idx]
        lam = lam.seq[lam.idx]
        x   = sim_pois_ar(n = n, phi = phi, lam = lam)

        # Estimate Poisson mean param
        lam.est = mean(x)

        # Calculate Hermite coefficients
        g.coefs = g_coefs(lam = lam.est, k = 1:n.coefs)

        # Coefficients of f: gam.x --> gam.z
        f.coefs = factorial(1:n.coefs) * g.coefs^2 / lam.est

        # coefficients of f^-1: gam.z --> gam.x
        finv.coefs = reversion(f.coefs)

        # sample acf X
        rho.x = acf(x, lag.max = 1, plot = FALSE, type = "correlation")$acf
        rho.x = rho.x[-1]

        # reversion to Z
        rho.z = power_series(rho.x, finv.coefs)

        # store param estimates matching previous format
        param.est = c(lam.est, rho.z)

        simResults.poisAR1[idx, "estim.method"] = "impliedYW"
        simResults.poisAR1[idx, "n"] = n
        simResults.poisAR1[idx, "lam.true"] = lam
        simResults.poisAR1[idx, "lam.est"]  = param.est[1]
        simResults.poisAR1[idx, "phi.true"] = phi
        simResults.poisAR1[idx, "phi.est"]  = param.est[2]
        idx = idx + 1
        setTxtProgressBar(pb, idx)
      }}}}
close(pb)

# Save results as .Rdata object
# simResults.poisAR1 <- simResults_poisson_ar1
# setwd("~/Desktop/LatentGaussCounts/tests/simulations/poisson-ar1")
# save(simResults_poisson_ar1, file = "simResults_poisson_ar1.Rdata")



