
# cdf for genpoisson
pgpois= function(x,lambda, eta) { # using Yisu's code for the cdf function
  cdf.vec <- rep(-99,length(x))
  for (i in 1:length(x)){
    cdf.vec[i] <- sum(dgpois(0:x[i],lambda,eta))
  }
  return(cdf.vec)
}


