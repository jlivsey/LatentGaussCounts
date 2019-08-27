#' Title
#'
#' @param theta
#' @param data
#'
#' @return
#' @export
#'
#'
likSIS_AR1 = function(theta, data, setseed=FALSE){

  if (setseed){
    set.seed(1)
  }
  cdf = function(x, lam){ ppois(q=x, lambda=lam) }
  pdf = function(x, lam){ dpois(x, lambda=lam) }

  theta1.idx = 1


  if (abs(theta[2])>=1){
    out = Inf
  }else{
    z.rest = function(a,b){
      # Generates N(0,1) variables restricted to (ai,bi),i=1,...n
      qnorm(runif(length(a),0,1)*(pnorm(b,0,1)-pnorm(a,0,1))+pnorm(a,0,1),0,1)
    }

    theta1 = theta[theta1.idx]
    n.theta1.idx = theta1.idx[length(theta1.idx)] # num params in theta1
    theta2.idx = (n.theta1.idx + 1):(n.theta1.idx + 1)
    phi = theta[theta2.idx]
    xt = data
    T1 = length(xt)
    N = 1000 # number of particles
    prt = matrix(0,N,T1) # to collect all particles
    wgh = matrix(0,N,T1) # to collect all particle weights

    a = qnorm(cdf(xt[1]-1,theta1),0,1)
    b = qnorm(cdf(xt[1],theta1),0,1)
    a = rep(a,N)
    b = rep(b,N)
    zprev = z.rest(a,b)
    zhat = phi*zprev
    prt[,1] = zhat

    wprev = rep(1,N)
    wgh[,1] = wprev

    for (t in 2:T1)
    {
      rt = sqrt(1-phi^2)
      a = (qnorm(cdf(xt[t]-1,theta1),0,1) - phi*zprev)/rt
      b = (qnorm(cdf(xt[t],theta1),0,1) - phi*zprev)/rt
      err = z.rest(a,b)
      znew = phi*zprev + rt*err
      zhat = phi*znew
      prt[,t] = zhat
      zprev = znew

      wgh[,t] = wprev*(pnorm(b,0,1) - pnorm(a,0,1))
      wprev = wgh[,t]
    }

    lik = pdf(xt[1],theta1)*mean(wgh[,T1])
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
