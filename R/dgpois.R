# pdf of generalized poisson
dgpois <- function(x,lambda,eta){
  tmp = log(lambda) + (x-1)*log(lambda+eta*x) -lgamma(x+1) - lambda - eta*x
  return(exp(tmp))
}
# # GenPoisson density from VGAM package--I need to check if this is right
# # code below uses lambda instead of eta and theta instead of eta

