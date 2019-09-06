#' Title
#'
#' @param theta
#' @param data
#'
#' @return
#' @export
#'
#'
likSIS_AR2 = function(theta, data, setseed=FALSE){

  if (setseed){
    set.seed(1)
  }
  cdf = function(x, lam){ ppois(q=x, lambda=lam) }
  pdf = function(x, lam){ dpois(x, lambda=lam) }

  theta1.idx = 1
  theta1 = theta[theta1.idx]
  n.theta1.idx = theta1.idx[length(theta1.idx)] # num params in theta1
  theta2.idx = (n.theta1.idx + 1):(n.theta1.idx + 2)
  phi = theta[theta2.idx]


  if (!((phi[2]+phi[1]<1) & (phi[2]-phi[1]<1) & (abs(phi[2])<1))){
    out = Inf
  }else{
    z.rest = function(a,b){
      # Generates N(0,1) variables restricted to (ai,bi),i=1,...n
      qnorm(runif(length(a),0,1)*(pnorm(b,0,1)-pnorm(a,0,1))+pnorm(a,0,1),0,1)
    }



    xt = data
    T1 = length(xt)
    N = 1000 # number of particles
    prt = matrix(0,N,T1) # to collect all particles
    wgh = matrix(0,N,T1) # to collect all particle weights

    a1 = qnorm(cdf(xt[1]-1,theta1),0,1)
    b1 = qnorm(cdf(xt[1],theta1),0,1)
    a1 = rep(a1,N)
    b1 = rep(b1,N)
    zprev1 = z.rest(a1,b1)
    zhat1 = phi[1]*zprev1
    prt[,1] = zhat1

    wprev = rep(1,N)
    wgh[,1] = wprev

    phi0 = if(abs(phi[1])<0.9){phi[1]}else{0.9*sign(phi[1])}
    rt2 = sqrt(1-phi0^2)
    a2 = (qnorm(cdf(xt[2]-1,theta1),0,1) - phi[1]*zprev1)/rt2
    b2 = (qnorm(cdf(xt[2],theta1),0,1) - phi[1]*zprev1)/rt2
    err2 = z.rest(a2,b2)
    znew2 = phi0*zprev1 + rt2*err2
    znew = c(znew2,zprev1)
    zhat2 = sum(phi*znew)
    prt[,2] = zhat2

    wgh[,2] = wprev*(pnorm(b2,0,1) - pnorm(a2,0,1))
    wprev = wgh[,2]

    zprev = znew

    for (t in 3:T1)
    {
      rt = 1/sqrt( ((1-phi[2])/(1+phi[2])) / ((1-phi[2])^2-phi[1]^2) )
      a = (qnorm(cdf(xt[t]-1,theta1),0,1) - sum(phi*zprev))/rt
      b = (qnorm(cdf(xt[t],theta1),0,1) - sum(phi*zprev))/rt
      err = z.rest(a,b)
      znew = sum(phi*zprev) + rt*err
      znew = c(znew,zprev[1])
      zhat = sum(phi*znew)
      prt[,t] = zhat
      zprev = znew

      wgh[,t] = wprev*(pnorm(b,0,1) - pnorm(a,0,1))
      wprev = wgh[,t]
    }

    lik = pdf(xt[1],theta1)*mean(na.omit(wgh[,T1]))
    nloglik = (-2)*log(lik)

    # out = list()
    # # out$lik = (if (is.na(nloglik)) Inf else nloglik)
    # out$lik = -lik
    # out$prt = prt
    # out$wgh = wgh

    # out =  (if (is.na(lik)) 0 else -lik)

    out = (if (is.na(nloglik) | lik==0) Inf else nloglik)
  }

  return(out)
}
