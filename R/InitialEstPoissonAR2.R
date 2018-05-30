# test example for initial estimation when Poisson and AR = 2
# I will generate some counts using our model with Poisson and AR = 2.
# Then I will fit AR(2) to the followingthree series: 1. The true Z's that
# I used to simulate the counts, 2. The count data, 3. The Z's I obtain from applying
# inverse transform to the counts. At first glance, and as probably expected
# it seems that the estimates from inverting the transform are ok especially when phi>0.

# simulation scheme
nsim = 1000
n = 100
lam = 2
# AR parameters
phi1 = 0.7
phi2 = -0.56
phi = c(phi1,phi2)

lhat = rep(0,nsim)
arTrue = matrix(,nsim,2)
ar = matrix(,nsim,2)
arCount = matrix(,nsim,2)
for (i in 1:nsim){
  #set.seed(Sys.time())
  ztrue = arima.sim(model = list(ar=phi), n = n); ztrue = ztrue/sd(ztrue) # standardized
  x = qpois(pnorm(ztrue), lam)

  # method of moments and reverse transform
  lhat[i] = mean(x)
  tmp = ppois(x,lhat[i])
  z = qnorm(tmp) # this doesnt take values on an interval

  arTrue[i,] = arima(ztrue,order = c(2,0,0),include.mean=0)$coef
  ar[i,] = arima(z,order = c(2,0,0),include.mean=0)$coef
  arCount[i,] = arima(x,order = c(2,0,0),include.mean=0)$coef
}

est = matrix(c(arTrue,ar,arCount),nsim,6)
boxplot(est, col=(c("red","red","darkgreen","darkgreen","blue","blue")),axes = F)
abline(v=2.5,col="gray",lty=2)
abline(v=4.5,col="gray",lty=2)
axis(1, at=1:6, labels=c("phi1", "phi2","phi1", "phi2","phi1", "phi2" ))
axis(2) #default way
box()
legend("topright", inset=.01,
       c("True","Inverse","Counts"), fill=c("red","darkgreen","blue"), horiz=F)



