
# ---- g() from equation (23) on 1-31-18 ----


#' g() function from latentGaussCounts paper. Equation (23) from draft
#' dated 1-28-18
#'
#' @param lam Poisson parameter
#' @param k function input parameter
#' @param polys list object containing the hermite polynomials. By default the
#' orthopolynorm package is used and preloaded Polys.Rdata is loaded with the
#' package
#'
#' @return function value
#'



g <- function(lam, k, polys=Polys){
  her <- as.function(polys[[k]]) # polys[[k]] = H_{k-1}
  N = which(round(ppois(1:100, lam), 7) == 1)[1]
  terms = exp(-qnorm(ppois(0:N, lam, lower.tail= TRUE))^2/2) *
    her(qnorm(ppois(0:N, lam, lower.tail = TRUE)))
  return(list(val=sum(terms)/sqrt(2*pi)/factorial(k), terms=terms))
}
