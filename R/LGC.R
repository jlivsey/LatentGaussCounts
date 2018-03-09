LGC <- function(x, count.family = "Poisson",
                   gauss.series = "AR",
                   estim.method = "gaussianLik",
                   max.terms = 30, p=NULL, d=NULL, q=NULL, n.mix=NULL, ...)
{
  # DEFINE a cdf function from count.family input. Make sure it accepts a vector
  #        valued input.
   if(count.family=="Poisson"){
     cdf = function(x, lam){ ppois(q=x, lambda=lam) }
     count.mean = function(lam){ lam }
     count.initial = function(data){ mean(data) }
     theta1.idx = 1
   }  else if(count.family=="mixed-Poisson"){
    if(is.null(n.mix)) stop("you must specify the number of Poissons to mix,
                             n.mix, to use count.family=mixed-Poisson")
    if(n.mix==1) stop("n.mix must be greater than 1")
    if(n.mix==2){
       cdf = function(x, theta){
         theta[1]*ppois(x, theta[2]) + (1-theta[1])*ppois(x, theta[3])
       }
       count.mean = function(theta){

       }
       count.initial = function(data){

       }
       theta1.idx = 1:(n.mix+1)
    }
    if(n.mix>2) stop("mixed-Poisson is not currently coded for n.mix>2")

  }else{ stop("please specify a valid count.family") }

  # DEFINE a gamz function from gauss.series input. Make sure it accepts a
  #        vector valued input.
  if(gauss.series=="AR"){
    if(is.null(p)) stop("you must specify the AR order, p, to use gauss.series=AR")
    if(p>1) stop("the ACVF is not coded for AR models of order higher than 1 currently")
    if(p==1){
      gamZ = function(h, phi){ phi^h / (1-phi^2) }
      gauss.initial = function(x){ acf(x, plot = FALSE)$acf[2] }
      theta2.idx = (theta1.idx+1):(theta1.idx+1)
    }else{ stop("the p specified is not valid") }
  }

  # ----------------------------------------------------------------------------

  if(estim.method=="gaussianLik"){
    g <- function(k, theta1, polys=Polys){
      her <- as.function(polys[[k]]) # polys[[k]] = H_{k-1}
      N = which(round(cdf(1:100, theta1), 7) == 1)[1]
      terms = exp(-qnorm(cdf(0:N, theta1))^2/2) *  her(qnorm(cdf(0:N, theta1)))
      return(sum(terms)/sqrt(2*pi)/factorial(k))
    }

    gamX = function(h, theta2, gamZ, g.vec, max.terms=30){
      sum(g.vec^2 * factorial(1:max.terms) * (gamZ(h, theta2))^(1:max.terms))
    }

    lik = function(theta, data){
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
      out = -2*dmvnorm(as.numeric(data), mean = mean.vec, sigma = Sigma, log = TRUE)
      return(out)
    }
  }

  initial.param = c(count.initial(x), gauss.initial(x))
  cat("initial parameter estimates: ", initial.param, "\n")
  optim.output <- optim(par = initial.param, fn = lik, data=x, ...)

  return(optim.output)
}
