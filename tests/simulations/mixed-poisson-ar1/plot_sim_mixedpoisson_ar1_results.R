load("simResults_mixedpoisAR1.Rdata")

D = simResults.mixedpoisAR1

op = par(mfrow=c(2,2), mar=c(2,2,4,1))

# Lambda estimates for positive phi values
d1 = density(D$lam1.est[D$n==200 & D$lam1.true==2 & D$phi.true==.75 & D$lam2.true==3])
d2 = density(D$lam2.est[D$n==200 & D$lam1.true==2 & D$phi.true==.75 & D$lam2.true==3])
plot(d2, xlim=c(0,7), main="lambda estimates (true=2)")
lines(d1, col=2)
abline(v=2, lty="dotted")
legend("topright", legend = c("400 - .75", "200 - .75", "100 - .75", "400 - .25", "200 - .25", "100 - .25"), col = c(1,2,3,1,2,3), lwd=c(1,1,1,3,3,3))


# Lambda estimates for positive phi values
d1 = density(D$lam1.est[D$n==200 & D$lam1.true==2 & D$phi.true==.75 & D$lam2.true==10])
d2 = density(D$lam2.est[D$n==200 & D$lam1.true==2 & D$phi.true==.75 & D$lam2.true==10])
plot(d2, xlim=c(0,20), main="lambda estimates (true=2)")
lines(d1, col=2)
abline(v=c(2, 10), lty="dotted")
legend("topright", legend = c("400 - .75", "200 - .75", "100 - .75", "400 - .25", "200 - .25", "100 - .25"), col = c(1,2,3,1,2,3), lwd=c(1,1,1,3,3,3))

par(op)
