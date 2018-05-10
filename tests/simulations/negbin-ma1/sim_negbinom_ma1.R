# setup params to sim from
n.seq   = c(200, 400)
tht.seq = c(-.75, -.25, .25, .75)
p.seq = c(.2, .5)
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
    x   = sim_negbin_ma1(n = n, tht=tht, p=p, r=r)
    out.optim = LGC(x,
                    count.family = "negbinom",
                    gauss.series = "MA", q=1,
                    estim.method = "gaussianLik")
    # store output in data.frame
    d$estim.method = "gaussianLik"
    d$n = n
    d$tht.true = tht
    d$tht.est  = out.optim$par[]
    d$tht.se   = out.ouptim$stderr[]
    d$p.true = p
    d$p.est = out.optim$par[]
    d$p.se = out.optim$stderr[]
    d$r.true = r
    d$r.est = out.optim$par[]
    d$r.se  = out.optim$stderr[]
    idx = idx + 1
    setTxtProgressBar(pb,idx)
}}}}}
close(pb)

simResults_negbin_ma1 = d
# save(simResults_negbin_ma1, file = "simResults_negbin_ma1.Rdata")