

x=sim_pois_ar(n = 400, phi = 0.7, lam = 2)
plot.ts(x)


lik_all = NULL
lik_all2 = NULL
phi_all = seq(.05,.95,.01)

for (phi in phi_all){
  lik = likSIS_AR1(c(2,phi),x)
  lik_all = c(lik_all,lik)
}

for (phi in phi_all){
  lik2 = likSIS_AR1(c(2,phi),x,TRUE)
  lik_all2 = c(lik_all2,lik2)
}

plot(phi_all,lik_all,type = "l")
lines(phi_all,lik_all2,type = "l",col="red")




count.mean = function(lam){ lam }
count.initial = function(data){ mean(data) }

gauss.initial = function(data){
    p=1
    lhat = mean(data)
    tmp = ppois(data,lhat)
    z = qnorm(tmp) # this doesnt take values on an interval but may still be used for initial AR estimates
    return(arima(z,order = c(p,0,0),include.mean=0)$coef)
}

initial.param = c(count.initial(x), gauss.initial(x))

optim.output <- optim(par = initial.param, fn = likSIS_AR1, setseed=TRUE,
                      data=x, hessian=TRUE, method = "BFGS")
optim.output




R = 100
phi.est.all = rep(0,R)
lam.est.all = rep(0,R)
phi.std.all = rep(0,R)
lam.std.all = rep(0,R)


for (r in 1:R){

  print(r)
  set.seed(r)
  x=sim_pois_ar(n = 400, phi = 0.7, lam = 2)
  initial.param = c(count.initial(x), gauss.initial(x))

  optim.output <- optim(par = initial.param, fn = likSIS_AR1, setseed=TRUE,
                        data=x, hessian=TRUE, method = "BFGS")

  phi.est.all[r]=optim.output$par[2]
  lam.est.all[r]=optim.output$par[1]
  phi.std.all[r]=sqrt(solve(optim.output$hessian))[2,2]
  lam.std.all[r]=sqrt(solve(optim.output$hessian))[1,1]

}


phi.est.up = phi.est.all + 1.96*phi.std.all*sqrt(2) # sqrt(2) since -2*loglik is minimized
phi.est.lo = phi.est.all - 1.96*phi.std.all*sqrt(2)
plot(seq(1,R),phi.est.all,ylim=c(min(phi.est.lo),max(phi.est.up)))
abline(h = .7, lty = 2, lwd = 1, col="red")
segments(1:R, phi.est.lo, 1:R, phi.est.up, col = "blue", lty = 1)


lam.est.up = lam.est.all + 1.96*lam.std.all*sqrt(2)
lam.est.lo = lam.est.all - 1.96*lam.std.all*sqrt(2)
plot(seq(1,R),lam.est.all,ylim=c(min(lam.est.lo),max(lam.est.up)))
abline(h = 2, lty = 2, lwd = 1, col="red")
segments(1:R, lam.est.lo, 1:R, lam.est.up, col = "blue", lty = 1)



###


xi1 = 2
xi2 = 5

xi1 = 2*(1+1i*sqrt(3))/3
xi2 = 2*(1-1i*sqrt(3))/3

phi1 = Re(1/xi1 + 1/xi2)
phi2 = -Re((1/xi1) * (1/xi2))

c(phi1,phi2)



x=sim_pois_ar(n = 200, phi = c(phi1,phi2), lam = 2)
plot.ts(x)
#acf(x,  lag.max = 20)

likSIS_ARp(c(2,phi1,phi2),x)
likSIS_AR2(c(2,phi1,phi2),x,setseed=TRUE)


res = likSIS_ARp(c(2,phi1,phi2),x)

count.mean = function(lam){ lam }
count.initial = function(data){ mean(data) }

gauss.initial = function(data){
  p=2
  lhat = mean(data)
  tmp = ppois(data,lhat)
  z = qnorm(tmp) # this doesnt take values on an interval but may still be used for initial AR estimates
  return(arima(z,order = c(p,0,0),include.mean=0)$coef)
}

initial.param = c(count.initial(x), gauss.initial(x))

optim.output <- optim(par = initial.param, fn = likSIS_ARp,
                      data=x, hessian=TRUE, method = "BFGS")
optim.output
