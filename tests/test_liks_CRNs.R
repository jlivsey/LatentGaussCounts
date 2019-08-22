

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


theta1.min = 0.01
theta1.max = mean(x) + 30

theta2.min = -.99
theta2.max = .99

optim.output2 <- optim(par = initial.param, fn = likSIS_AR1,
                      data=x, setseed=TRUE, hessian=TRUE, method = "L-BFGS",
                      lower = c(theta1.min, theta2.min),
                      upper = c(theta1.max, theta2.max))
optim.output2





