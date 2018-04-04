load("simResults_poisAR1.Rdata")

D = data.frame(n = as.numeric(simResults.poisAR1[,"n"]),
               lam.true = as.numeric(simResults.poisAR1[, "lam.true"]),
               lam.est = as.numeric(simResults.poisAR1[, "lam.est"]),
               phi.true = as.numeric(simResults.poisAR1[, "phi.true"]),
               phi.est = as.numeric(simResults.poisAR1[, "phi.est"]))


op = par(mfrow=c(2,2), mar=c(2,2,4,1))

# Lambda estimates for positive phi values
d1 = density(D$lam.est[A$n==100 & D$lam.true==2 & D$phi.true==.75])
d2 = density(D$lam.est[A$n==200 & D$lam.true==2 & D$phi.true==.75])
d4 = density(D$lam.est[A$n==400 & D$lam.true==2 & D$phi.true==.75])
d1s = density(D$lam.est[A$n==100 & D$lam.true==2 & D$phi.true==.25])
d2s = density(D$lam.est[A$n==200 & D$lam.true==2 & D$phi.true==.25])
d4s = density(D$lam.est[A$n==400 & D$lam.true==2 & D$phi.true==.25])
plot(d4s, xlim=c(1, 3.5), main="lambda estimates (true=2)")
lines(d4, col=1)
lines(d2, col=2)
lines(d1, col=3)
lines(d4s, col=1, lwd=3)
lines(d2s, col=2, lwd=3)
lines(d1s, col=3, lwd=3)
abline(v=2, lty="dotted")
legend("topright", legend = c("400 - .75", "200 - .75", "100 - .75", "400 - .25", "200 - .25", "100 - .25"), col = c(1,2,3,1,2,3), lwd=c(1,1,1,3,3,3))

# Lambda estimates for negative phi values
d1 = density(D$lam.est[A$n==100 & D$lam.true==2 & D$phi.true==-.75])
d2 = density(D$lam.est[A$n==200 & D$lam.true==2 & D$phi.true==-.75])
d4 = density(D$lam.est[A$n==400 & D$lam.true==2 & D$phi.true==-.75])
d1s = density(D$lam.est[A$n==100 & D$lam.true==2 & D$phi.true==-.25])
d2s = density(D$lam.est[A$n==200 & D$lam.true==2 & D$phi.true==-.25])
d4s = density(D$lam.est[A$n==400 & D$lam.true==2 & D$phi.true==-.25])
plot(d4, xlim=c(1, 3.5), main="lambda estimates (true=2)")
lines(d4, col=1)
lines(d2, col=2)
lines(d1, col=3)
lines(d4s, col=1, lwd=3)
lines(d2s, col=2, lwd=3)
lines(d1s, col=3, lwd=3)
abline(v=2, lty="dotted")
legend("topright", legend = c("400 - -.75", "200 - -.75", "100 - -.75", "400 - -.25", "200 - -.25", "100 - -.25"), col = c(1,2,3,1,2,3), lwd=c(1,1,1,3,3,3))

# positive phi estimates of phi.true=.75 for lambda=2 and 10
d1 = density(D$phi.est[A$n==100 & D$lam.true==10 & D$phi.true==.75])
d2 = density(D$phi.est[A$n==200 & D$lam.true==10 & D$phi.true==.75])
d4 = density(D$phi.est[A$n==400 & D$lam.true==10 & D$phi.true==.75])
d1s = density(D$phi.est[A$n==100 & D$lam.true==2 & D$phi.true==.75])
d2s = density(D$phi.est[A$n==200 & D$lam.true==2 & D$phi.true==.75])
d4s = density(D$phi.est[A$n==400 & D$lam.true==2 & D$phi.true==.75])
plot(d4s, main="Phi estimates")
lines(d4, col=1)
lines(d2, col=2)
lines(d1, col=3)
lines(d4s, col=1, lwd=3)
lines(d2s, col=2, lwd=3)
lines(d1s, col=3, lwd=3)
abline(v=.75, lty="dotted")
legend("topright", legend = c("400 - 10", "200 - 10", "100 - 10", "400 - 2", "200 - 2", "100 - 2"), col = c(1,2,3,1,2,3), lwd=c(1,1,1,3,3,3))

# positive phi estimates of phi.true=.25 for lambda=2 and 10
d1 = density(D$phi.est[A$n==100 & D$lam.true==10 & D$phi.true==-.75])
d2 = density(D$phi.est[A$n==200 & D$lam.true==10 & D$phi.true==-.75])
d4 = density(D$phi.est[A$n==400 & D$lam.true==10 & D$phi.true==-.75])
d1s = density(D$phi.est[A$n==100 & D$lam.true==2 & D$phi.true==-.75])
d2s = density(D$phi.est[A$n==200 & D$lam.true==2 & D$phi.true==-.75])
d4s = density(D$phi.est[A$n==400 & D$lam.true==2 & D$phi.true==-.75])
plot(d4, main="Phi estimates")
lines(d4, col=1)
lines(d2, col=2)
lines(d1, col=3)
lines(d4s, col=1, lwd=3)
lines(d2s, col=2, lwd=3)
lines(d1s, col=3, lwd=3)
abline(v=-.75, lty="dotted")
legend("topright", legend = c("400 - 10", "200 - 10", "100 - 10", "400 - 2", "200 - 2", "100 - 2"), col = c(1,2,3,1,2,3), lwd=c(1,1,1,3,3,3))

par(op)
