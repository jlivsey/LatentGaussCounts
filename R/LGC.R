#' Latent Gaussian Count model builder
#'
#' Main execution function for latentGaussCounts package.
#'
#' @param x data
#' @param count.family desired marginal distribution (Poisson,mixed-Poisson,etc)
#' @param gauss.series desired structure of your latent Gaussian process
#' @param estim.method method used for estimating parameters
#' @param max.terms maximum number of terms used to truncate Hermite expansions
#' @param p AR order
#' @param d don't use it
#' @param q don't use it
#' @param n.mix number of Poisson distributions in mixed-Poisson count.family
#' @param n Not sure what this is
#' @param print.progress Should progress be printed as optim() is run
#' @param print.initial.estimates Should initial estiamtes be printed
#' @param ... additional parameters to pass to optim
#'
#' @return list object which are results from optim() run
#'
#' @export

LGC <- function(x, count.family = c("Poisson", "mixed-Poisson", "negbinom", "GenPoisson"),
                   gauss.series = c("AR","FARIMA"),
                   estim.method = c("gaussianLik","particlesSIS"),
                   max.terms = 30, p=NULL, d=NULL, q=NULL, n.mix=NULL, n=NULL,
                   print.progress=FALSE, print.initial.estimates=FALSE, ...)
{

  # ---- count.family input ---------------------------------------------------

  # DEFINE a cdf function from count.family input. Make sure it accepts a vector
  #        valued input.
   if(count.family=="Poisson"){
     cdf = function(x, lam){ ppois(q=x, lambda=lam) }
     pdf = function(x, lam){ dpois(x, lambda=lam) }
     count.mean = function(lam){ lam }
     count.initial = function(data){ mean(data) }
     theta1.min = 0.01
     theta1.max = mean(x) + 30
     theta1.idx = 1
   }
  #----------------------------------------------------------------------------------------------#
  #-------------------------------------Stef--May 22---------------------------------------------#
    else if(count.family=="GenPoisson"){
     cdf= function(x,parameter) {
       lambda = parameter[1]
       theta = parameter[2]
       cdf.vec <- rep(-99,length(x))
       for (i in 1:length(x)){
         cdf.vec[i] <- sum(pdf(0:x[i],parameter))
       }
       return(cdf.vec)
     }
     # I am using the source code of the dgenpois function from the vglm package
     pdf = function(x, parameter){
       lambda = parameter[1]
       theta = parameter[2]
       dgenpois(x, lambda, theta, log = FALSE)
     }
     count.mean = function(parameter){
       lambda = parameter[1]
       theta = parameter[2]
       return(theta/(1-lambda))
     }
     # FIX ME: I ll use method of moments to get initial values but this will only work when
     # lambda i between 0 and 1
     count.initial = function(data){
       lambda.hat = 1 - sqrt(mean(data))/sd(data) # FIX ME: solving for variance yields tqo solution plus/minus
                                             # need to compute likelihood for both and select the best
       theta.hat = mean(data)*(1-lambda.hat)
       return(c(lambda.hat, theta.hat))
     }
     theta1.min = c(0.01,0.01)
     theta1.max = c(0.99,mean(x) + 30)
     theta1.idx = 2 # FIX ME: I am not sure what is this
    #----------------------------------------------------------------------------------------------#
   }else if(count.family=="mixed-Poisson"){
    if(is.null(n.mix)) stop("you must specify the number of Poissons to mix,
                             n.mix, to use count.family=mixed-Poisson")
    if(n.mix==1) stop("n.mix must be greater than 1")
    if(n.mix==2){
       cdf = function(x, theta){
         theta[1]*ppois(x, theta[2]) + (1-theta[1])*ppois(x, theta[3])
       }
       pdf = function(x, theta){
         theta[1]*dpois(x, theta[2]) + (1-theta[1])*dpois(x, theta[3])
       }
       count.mean = function(theta){
          theta[1]*theta[2] + (1-theta[1])*theta[3]
       }
       count.initial = function(data){
         n = length(data)
         n2 = floor(length(data)/2)
         p.hat = 1/2
         theta1.hat = mean(sort(data)[1:n2])
         theta2.hat = mean(sort(data)[(n2+1):n])
         return(c(p.hat, theta1.hat, theta2.hat))
       }
       theta1.min = c(0,0.01,0.01)
       theta1.max = c(0.5,20,20)
       theta1.idx = 1:(n.mix+1)
    }
    if(n.mix>2) stop("mixed-Poisson is not currently coded for n.mix>2")

   } else if(count.family=="Binomial"){
     if(is.null(n)) stop("you must specify the number of trials,
                             n, to use count.family=Binomial")
     cdf = function(x, theta){ pbinom(q = x, size = n, prob = theta) }
     pdf = function(x, theta){ dbinom(x = x, size = n, prob = theta) }
     count.mean = function(theta){ n*theta }
     count.initial = function(data){ mean(data)/n }
     theta1.idx = 1
   } else if(count.family=="negbinom"){
     cdf = function(x, theta){
       pnbinom(q = x, size = theta[1], prob = theta[2])
     }
     pdf = function(x, theta){
       dnbinom(x = x, size = theta[1], prob = theta[2])
     }
     count.mean = function(theta){
       theta[1] * (1-theta[2]) / theta[2]
     }
     count.initial = function(data){
       n = length(data)
       m = mean(x)
       v = var(x)
       theta1.hat = m^2 / (v-m)
       theta2.hat = m/v
       return(c(theta1.hat, theta2.hat))
     }
     theta1.min = c(0.01, 0.1)
     theta1.max = c(10  , 0.9)
     theta1.idx = 1:2
   }else{ stop("please specify a valid count.family") }


  # ---- Gauss.series input ---------------------------------------------------

  # DEFINE a gamz function from gauss.series input. Make sure it accepts a
  #        vector valued input.
  if(gauss.series=="AR"){
    if(is.null(p)) stop("you must specify the AR order, p, to use
                        gauss.series=AR")
    if(p>3) stop("the ACVF is not coded for AR models of order higher than 1
                 currently")
    if(p==1){
      gamZ = function(h, phi){ phi^h}
      gauss.initial = function(x){ acf(x, plot = FALSE)$acf[2] }
      n.theta1.idx = theta1.idx[length(theta1.idx)] # num params in theta1
      theta2.idx = (n.theta1.idx + 1):(n.theta1.idx + 1)
      theta2.min = -.99
      theta2.max = .99
    } else if(p==2){
      gamZ = function(h, phi){ ARMAacf(ar = phi, lag.max = 1000)[h+1] }
      gauss.initial = function(data){ acf(data, plot = FALSE)$acf[2:3] }
      n.theta1.idx = theta1.idx[length(theta1.idx)] # num params in theta1
      theta2.idx = (n.theta1.idx + 1):(n.theta1.idx + 2)
    } else if(p==3){
      gamZ = function(h, phi){ ARMAacf(ar = phi, lag.max = 1000)[h+1] }
      gauss.initial = function(data){ acf(data, plot = FALSE)$acf[2:4] }
      n.theta1.idx = theta1.idx[length(theta1.idx)] # num params in theta1
      theta2.idx = (n.theta1.idx + 1):(n.theta1.idx + 3)
    } else{ stop("the p specified is not valid") }
  }

  if(gauss.series=="FARIMA"){ # currently only FARIMA(0,d,0)
    gamZ = function(h, d){ acf.farima0d0(d = d, h = h) }
    gauss.initial = function(data){ 1/4 }
    n.theta1.idx = theta1.idx[length(theta1.idx)] # num params in theta1
    theta2.idx = (n.theta1.idx + 1):(n.theta1.idx + 1)
  }

  if(gauss.series=="MA"){
    gamZ = function(h, tht){
        if(h==0){
            return(1)
      } else if(h==1){
            return(tht/(1+tht^2))
      } else{
            return(0)
      }
    }
    gauss.initial = function(x){ acf(x, plot = FALSE)$acf[2] }
    n.theta1.idx = theta1.idx[length(theta1.idx)] # num params in theta1
    theta2.idx = (n.theta1.idx + 1):(n.theta1.idx + 1)
    theta2.min = -.99
    theta2.max = .99
  }

# ---- Estimation Methods -----------------------------------------------------

  if(estim.method=="gaussianLik"){
    g <- function(k, theta1){
      #her <- as.function(Polys[[k]]) # polys[[k]] = H_{k-1}
      N = which(round(cdf(1:10000, theta1), 7) == 1)[1]
      if(length(N)==1 | is.na(N) ) stop("Haven't reached upper limit for cdf")
      terms = exp(-qnorm(cdf(0:N, theta1))^2/2) *
              Hermite_poly(k = k, x = qnorm(cdf(0:N, theta1)))
             #her(qnorm(cdf(0:N, theta1)))
      return(sum(terms)/sqrt(2*pi)/factorial(k))
    }

    gamX = function(h, theta2, gamZ, g.vec, max.terms=30){
      sum(g.vec^2 * factorial(1:max.terms) * (gamZ(h, theta2))^(1:max.terms))
    }

    lik = function(theta, data){
      if(print.progress) cat("theta = ", theta, "\n")
      theta1 = theta[theta1.idx]
      theta2 = theta[theta2.idx]
      n   = length(data)
      h   = 0:(n-1)
      g.vec = c()
      for(k in 1:30){g.vec[k] <- g(k=k, theta1 = theta1)}
      gamX.vec = c()
      for(i in 1:length(h)){gamX.vec[i] = gamX(h[i], theta2, gamZ, g.vec)}
      Sigma = toeplitz(gamX.vec)
      mean.vec = rep(count.mean(theta1), n)
      out = -2*mvtnorm::dmvnorm(as.numeric(data), mean = mean.vec,
                                sigma = Sigma, log = TRUE)
      if(print.progress) cat(" lik = ", out, "\n")
      if(out=="Inf"){ return(10^6) }
      return(out)
    }
  }

  if ((gauss.series=="AR") & (estim.method=="particlesSIS")){
    if(is.null(p)) stop("you must specify the AR order, p, to use
                        gauss.series=AR")
    if(p>2) stop("the ACVF is not coded for AR models of order higher than 1
                 currently")
    if(p==1){
      #set.seed(1)
      z.rest = function(a,b){
        # Generates N(0,1) variables restricted to (ai,bi),i=1,...n
        qnorm(runif(length(a),0,1)*(pnorm(b,0,1)-pnorm(a,0,1))+pnorm(a,0,1),0,1)
      }
      likSIS = function(theta, data){
        #set.seed(1)
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

        out = if (is.na(nloglik)) Inf else nloglik
        return(out)
      }
    } else if(p==2){
      #set.seed(1)
      z.rest = function(a,b){
        # Generates N(0,1) variables restricted to (ai,bi),i=1,...n
        qnorm(runif(length(a),0,1)*(pnorm(b,0,1)-pnorm(a,0,1))+pnorm(a,0,1),0,1)
      }
      likSIS = function(theta, data){
        # set.seed(1)
        theta1 = theta[theta1.idx]
        n.theta1.idx = theta1.idx[length(theta1.idx)] # num params in theta1
        theta2.idx = (n.theta1.idx + 1):(n.theta1.idx + 2)
        phi = theta[theta2.idx]

        if ( (phi[2]+phi[1]<1) & (phi[2]-phi[1]<1) & (abs(phi[2])<1)){
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

          lik = pdf(xt[1],theta1)*mean(wgh[,T1])
          nloglik = (-2)*log(lik)

          out = if (is.na(nloglik)) Inf else nloglik
          return(out)
        }else{
          out = Inf
          return(out)
        }
      }
    } else{ stop("the p specified is not valid") }
  }

  if ((gauss.series=="FARIMA") & (estim.method=="particlesSIS")){
    #set.seed(1)
    z.rest = function(a,b){
      # Generates N(0,1) variables restricted to (ai,bi),i=1,...n
      qnorm(runif(length(a),0,1)*(pnorm(b,0,1)-pnorm(a,0,1))+pnorm(a,0,1),0,1)
    }
    likSIS = function(theta, data){
      #set.seed(1)
      theta1 = theta[theta1.idx]
      n.theta1.idx = theta1.idx[length(theta1.idx)] # num params in theta1
      theta2.idx = (n.theta1.idx + 1):(n.theta1.idx + 1)
      d = theta[theta2.idx]
      xt = data
      N = 1000 # number of particles

      T1 = length(xt)
      prt = matrix(0,N,T1) # to collect all particles
      wgh = matrix(0,N,T1) # to collect all particle weights

      a = qnorm(cdf(xt[1]-1,theta1),0,1)
      b = qnorm(cdf(xt[1],theta1),0,1)
      a = rep(a,N)
      b = rep(b,N)
      zprev = z.rest(a,b)
      phi1 = -(gamma(2)/gamma(2-d))*(gamma(1-d)/gamma(2))*(gamma(1-d)/gamma(1))*(if (d==0) {0} else {1/gamma(-d)})
      zhat = phi1*rev(zprev)
      prt[,1] = zhat

      wprev = rep(1,N)
      wgh[,1] = wprev

      rt = 1

      for (t in 2:T1)
      {
        rt = rt*sqrt(1-(d/(t-1-d))^2)
        a = (qnorm(cdf(xt[t]-1,theta1),0,1) - zhat)/rt
        b = (qnorm(cdf(xt[t],theta1),0,1) - zhat)/rt
        err = z.rest(a,b)
        znew = zhat + rt*err
        zprev = cbind(zprev,znew)

        phit = -exp(lgamma(t+1)-lgamma(t-d+1))*exp(lgamma((1:t)-d)-lgamma((1:t)+1))*exp(lgamma(t-d-(1:t)+1)-lgamma(t-(1:t)+1))*(if (d==0) {0} else {1/gamma(-d)})
        zhat = rowSums(phit*zprev[,ncol(zprev):1])
        prt[,t] = zhat

        wgh[,t] = wprev*(pnorm(b,0,1) - pnorm(a,0,1))
        wprev = wgh[,t]
      }

      lik <- dpois(xt[1],theta1)*mean(wgh[,T1])
      # lik = dpois(xt[1],theta1)*mean(wgh[,T1],na.rm = TRUE)
      nloglik = (-2)*log(lik)

      out = if (is.na(nloglik)) Inf else nloglik
      return(out)
    }
  }

  # ---- Optimization ---------------------------------------------------------

  if(estim.method=="gaussianLik"){
    initial.param = c(count.initial(x), gauss.initial(x))
    if(print.initial.estimates) cat("initial parameter estimates: ", initial.param, "\n")
    optim.output <- optim(par = initial.param, fn = lik,
                          data=x, hessian=TRUE, method = "L-BFGS-B",
                          lower = c(theta1.min, theta2.min),
                          upper = c(theta1.max, theta2.max))
    stder = rep(NA, length(initial.param))
    stder <- sqrt(diag(solve(optim.output$hessian))) # calculate standard error
    stder[is.nan(stder)] = 0
    optim.output = append(optim.output, list(stder=stder)) # append stder output
  }

  if((gauss.series=="AR") & (estim.method=="particlesSIS")){
    if(p==1){
      R0 <- DEoptim::DEoptim(likSIS, lower = c(theta1.min,-.99), upper = c(theta1.max,.99),control = DEoptim::DEoptim.control(trace = 10, itermax = 100, steptol = 50, reltol = 1e-5), data=x)
      optim.output <- as.vector(R0$optim$bestmem)
    }
    if(p==2){
      R0 <- DEoptim::DEoptim(likSIS, lower = c(theta1.min,-1.99,-.99), upper = c(theta1.max,1.99,.99),control = DEoptim::DEoptim.control(trace = 10, itermax = 100, steptol = 50, reltol = 1e-5), data=x)
      optim.output <- as.vector(R0$optim$bestmem)
    }
  }

  if((gauss.series=="FARIMA") & (estim.method=="particlesSIS")){
      R0 <- DEoptim::DEoptim(likSIS, lower = c(theta1.min,-.49), upper = c(theta1.max,.49),control = DEoptim::DEoptim.control(trace = 10, itermax = 100, steptol = 50, reltol = 1e-5), data=x)
      optim.output <- as.vector(R0$optim$bestmem)
  }

  return(optim.output)
}
