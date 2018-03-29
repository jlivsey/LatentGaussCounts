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
