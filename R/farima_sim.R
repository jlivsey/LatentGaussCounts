
#' Simulate from FARIMA(0,d,0)
#'
#' @param Tt length of the time series to be generated
#' @param d long memory parameter
#'
#' @return time series of length Tt
#' @export
#'
#' @examples farima.sim(200, 1/4)


farima.sim <- function(Tt,d){
  r <- rep(0,2*Tt-2)
  r[1] <- gamma(1-2*d)/(gamma(1-d)^2)
  r[2:Tt] <- gamma(1-2*d)/gamma(d)/gamma(1-d)*exp(lgamma((1:(Tt-1))+d)-lgamma((1:(Tt-1))-d+1))
  r[(Tt+1):(2*Tt-2)] <- rev(r[2:(Tt-1)])
  # James's embedding
  L <- Re(fft(r))
  wr <- rnorm(2*Tt-2,0,1)
  wi <- rnorm(2*Tt-2,0,1)
  w <- wr + 1i*wi
  w2 <- sqrt(L/(2*Tt-2))*w

  z <- fft(w2)
  return(Re(z[1:Tt]))
}
