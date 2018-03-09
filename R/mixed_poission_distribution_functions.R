# .........................................................
# Mixed Poisson distribution functions
# .........................................................


rmixpois = function(n, p, lam1, lam2){
  u = runif(n)
  x = rpois(n,lam1)*(u<=p) + rpois(n,lam2)*(u>p)
  return(x)
}

dmixpois = function(x, p, lam1, lam2){
  y = p*dpois(x,lam1) + (1-p)*dpois(x,lam2)
  return(y)
}

pmixpois = function(x, p, lam1, lam2){
  y = p*ppois(x,lam1) + (1-p)*ppois(x,lam2)
  return(y)
}

qmixpois = function(y, p, lam1, lam2){
  yl = length(y)
  x = rep(0,yl)
  for (n in 1:yl){
    while(pmixpois(x[n], p, lam1, lam2) <= y[n]){ # R qpois would use <y; this choice makes the function right-continuous; this does not really matter for our model
      x[n] = x[n]+1
    } 
  }
  return(x)
}


# .........................................................
# Generate mixed Poisson-AR1 model
# Inputs:
#  T: length of the series
#  p: mixture porbability
#  lam1, lam2: Poisson lambda parameters
#  phi: AR(1) phi parameter
# .........................................................

mixpois.ar1 <- function(Tt,phi,p,lam1,lam2){
  zt <- rep(0,Tt)
  zt[1] <- rnorm(1,0,1)
  for (t in 2:Tt){
    zt[t] = phi*zt[t-1] + sqrt(1-phi^2)*rnorm(1,0,1)
  }
  xt <- qmixpois(pnorm(zt,0,1),p,lam1,lam2)
  return(xt)
}