# setup params to sim from
n.seq   = c(200)
tht.seq = c(-.75, .75)
p.seq = c(.2)
r.seq = 3
Nsim = 200 # total iterations per param configuration

total.iter = Nsim * length(n.seq) * length(tht.seq) *
             length(p.seq) * length(r.seq)
save.cols = c("estim.method", "n",
              "tht.true", "tht.est", "tht.se",
              "p.true", "p.est", "p.se",
              "r.true", "r.est", "r.se")
d = data.frame(matrix(nrow = total.iter, ncol = length(save.cols)))
colnames(d) = save.cols

print(sprintf("total iterations: %s", total.iter))
idx = 1
pb <- txtProgressBar(min=2,max=total.iter,style=3)
for(n.idx in 1:length(n.seq)){
for(tht.idx in 1:length(tht.seq)){
for(p.idx in 1:length(p.seq)){
for(r.idx in 1:length(r.seq)){
for(dumb.variable in 1:Nsim){
    n   = n.seq[n.idx]
    tht = tht.seq[tht.idx]
    p   = p.seq[p.idx]
    r   = r.seq[r.idx]
    x   = sim_negbinom_ma1(n = n, theta=tht, p=p, r=r)
    out.optim = LGC(x,
                    count.family = "negbinom",
                    gauss.series = "MA", q=1,
                    estim.method = "gaussianLik")
    # store output in data.frame
    d$estim.method[idx] = "gaussianLik"
    d$n = n[idx]
    d$tht.true[idx] = tht
    d$tht.est[idx]  = out.optim$par[3]
    d$tht.se[idx]   = out.optim$stder[3]
    d$p.true[idx] = p
    d$p.est[idx] = out.optim$par[2]
    d$p.se[idx] = out.optim$stder[2]
    d$r.true[idx] = r
    d$r.est[idx] = out.optim$par[1]
    d$r.se[idx]  = out.optim$stder[1]
    idx = idx + 1
    setTxtProgressBar(pb,idx)
}}}}}
close(pb)

simResults_negbin_ma1 = d
# save(simResults_negbin_ma1, file = "simResults_negbin_ma1.Rdata")
