#' Title
#'
#' @param theta
#' @param data
#'
#' @return
#' @export
#'
#'
likSIS_ARp = function(theta, data){

  cdf = function(x, lam){ ppois(q=x, lambda=lam) }
  pdf = function(x, lam){ dpois(x, lambda=lam) }


  out = list()
  set.seed(1)

  p = length(theta)-1
  z.rest = function(a,b){
    # Generates N(0,1) variables restricted to (ai,bi),i=1,...n
    qnorm(runif(length(a),0,1)*(pnorm(b,0,1)-pnorm(a,0,1))+pnorm(a,0,1),0,1)
  }


  theta1.idx = 1

  theta1 = theta[theta1.idx]
  n.theta1.idx = theta1.idx[length(theta1.idx)] # num params in theta1
  theta2.idx = (n.theta1.idx + 1):(n.theta1.idx + p)
  phi = theta[theta2.idx]


  if (prod(abs(polyroot(c(1,-phi))) > 1)){ # that is the ar model is causal

    xt = data
    T1 = length(xt)
    N = 1000 # number of particles
    prt = matrix(0,N,T1) # to collect all particles
    wgh = matrix(0,N,T1) # to collect all particle weights

    library(FitAR) # this probably needs to be uploaded differently

    circulant <- function(x, nrow = length(x)) {
      n <- length(x)
      matrix(x[(1:n - rep(1:nrow, each=n)) %% n + 1L], ncol=n, byrow=TRUE)
    }

    a = qnorm(cdf(xt[1]-1,theta1),0,1)
    b = qnorm(cdf(xt[1],theta1),0,1)
    a = rep(a,N)
    b = rep(b,N)
    zprev = z.rest(a,b)
    zhat = (TacvfAR(phi)[2]/TacvfAR(phi)[1])*zprev
    prt[,1] = zhat

    wprev = rep(1,N)
    wgh[,1] = wprev

    if (p==1){
      rt = sqrt(1-phi^2)
    }

    if (p>=2){
      for (t in 2:p){

        Gt = circulant(TacvfAR(phi)[1:t])
        gt = TacvfAR(phi)[2:(t+1)]
        phit = solve(Gt) %*% gt
        rt =  sqrt(TacvfAR(phi)[1] - gt %*% solve(Gt) %*% gt)

        a = (qnorm(cdf(xt[t]-1,theta1),0,1) - sum(phit*zprev))/rt
        b = (qnorm(cdf(xt[t],theta1),0,1) - sum(phit*zprev))/rt
        err = z.rest(a,b)
        znew = sum(phti*zprev) + rt*err
        znew = c(znew,zprev[1])
        zhat = sum(phit*znew)
        prt[,t] = zhat
        zprev = znew

        wgh[,t] = wprev*(pnorm(b,0,1) - pnorm(a,0,1))
        wprev = wgh[,t]
      }
    }

    for (t in (p+1):T1){

      if (p==1){
        a = (qnorm(cdf(xt[t]-1,theta1),0,1) - phi*zprev)/rt
        b = (qnorm(cdf(xt[t],theta1),0,1) - phi*zprev)/rt
        err = z.rest(a,b)
        znew = phi*zprev + rt*err
        zhat = phi*znew
        prt[,t] = zhat
        zprev = znew
      }

      if (p>=2){
        a = (qnorm(cdf(xt[t]-1,theta1),0,1) - sum(phi*zprev))/rt
        b = (qnorm(cdf(xt[t],theta1),0,1) - sum(phi*zprev))/rt
        err = z.rest(a,b)
        znew = sum(phi*zprev) + rt*err
        if (p>=2){
          znew = c(znew,zprev[p-1])
        }
        zhat = sum(phi*znew)
        prt[,t] = zhat
        zprev = znew
      }

      wgh[,t] = wprev*(pnorm(b,0,1) - pnorm(a,0,1))
      wprev = wgh[,t]
    }

    lik = pdf(xt[1],theta1)*mean(wgh[,T1])
    nloglik = (-2)*log(lik)
    out$lik = (if (is.na(nloglik) | lik==0) Inf else nloglik)
    out$prt =prt
    out$wgh = wgh

  }else{
    out$lik = Inf
  }

  return(out)
}
