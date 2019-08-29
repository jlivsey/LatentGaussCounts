n.seq    = c(100, 200, 400)
phi.seq  = c(-.75, -.25, .25, .75)
p.seq    = 1/4
lam1.seq = c(2)
lam2.seq = c(5, 10)
Nsim = 200
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


# --- Add Gaussian Likelihood results to simResults matrix ----

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


# ---- Add Particle filtering results to SimResults matrix ----
# ---- Windows parallel -------------------------------------------------------
library(parallel)
# load data as list
l2 <- list()
for(i in 1:100) l2[[i]] = sim_mixedpois_ar1(n = 200, phi = .7, p = .25, lam1 = 2, lam2 = 10)
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
                                          count.family = "mixed-Poisson", n.mix=2,
                                          gauss.series = "AR", p=1,
                                          estim.method = "particlesSIS")
                                    })
})
# This is first 5 runs
system.time({
out2 <- parLapply(cl, l2[6:25], function(x){
                                      LGC(x,
                                          count.family = "mixed-Poisson", n.mix=2,
                                          gauss.series = "AR", p=1,
                                          estim.method = "particlesSIS")
                                    })
})
# This is next 20 runs (add previous 2 for first 25 runs)
system.time({
out3 <- parLapply(cl, l2[26:100], function(x){
                                      LGC(x,
                                          count.family = "mixed-Poisson", n.mix=2,
                                          gauss.series = "AR", p=1,
                                          estim.method = "particlesSIS")
                                    })
})
# Stop the cluster
stopCluster(cl)


# ---- Add implied Yule-Walker estim results to SimResults matrix ----
library(latentGaussCounts)
library(orthopolynom)
polys <- hermite.he.polynomials(100) # In global env for g_coefs() function - BAD PROGRAMING
source("~/Dropbox/jim/MVcopula/code/functions.R")
n.coefs <- 20

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

            # Estimate Mixed-Poisson params
            n = length(x)
            p.hat = 1/4
            n2 = floor(length(x)*p.hat)
            theta1.hat = mean(sort(x)[1:n2])
            theta2.hat = mean(sort(x)[(n2+1):n])
            param.est <- c(p.hat, theta1.hat, theta2.hat)

            # Calculate Hermite coefficients
            g.coefs = g_coefs_cdf(pmixpois_vecparam, param.est, k = 1:n.coefs)

            # Coefficients of f: gam.x --> gam.z
            process.var <- var(x)
            f.coefs = factorial(1:n.coefs) * g.coefs^2 / process.var

            # coefficients of f^-1: gam.z --> gam.x
            finv.coefs = reversion(f.coefs)

            # sample acf X
            rho.x = acf(x, lag.max = 1, plot = FALSE, type = "correlation")$acf
            rho.x = rho.x[-1]

            # reversion to Z
            rho.z = power_series(rho.x, finv.coefs)

            simResults.mixedpoisAR1$estim.method[idx] = "impliedYW"
            simResults.mixedpoisAR1$n[idx]         = n
            simResults.mixedpoisAR1$p.true[idx]    = p
            simResults.mixedpoisAR1$p.est[idx]     = p.hat
            simResults.mixedpoisAR1$p.se[idx]      = NA
            simResults.mixedpoisAR1$lam1.true[idx] = lam1
            simResults.mixedpoisAR1$lam1.est[idx]  = theta1.hat
            simResults.mixedpoisAR1$lam1.se[idx]   = NA
            simResults.mixedpoisAR1$lam2.true[idx] = lam2
            simResults.mixedpoisAR1$lam2.est[idx]  = theta2.hat
            simResults.mixedpoisAR1$lam2.se[idx]   = NA
            simResults.mixedpoisAR1$phi.true[idx]  = phi
            simResults.mixedpoisAR1$phi.est[idx]   = rho.z
            simResults.mixedpoisAR1$phi.se[idx]    = NA
            idx = idx + 1
            setTxtProgressBar(pb,idx)
          }}}}}}
close(pb)

# Save results as .Rdata object
# simResults.poisAR1 <- simResults_poisson_ar1
# setwd("~/Desktop/LatentGaussCounts/tests/simulations/poisson-ar1")
# save(simResults_poisson_ar1, file = "simResults_poisson_ar1.Rdata")



