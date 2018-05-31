PDFdgenpois <- function(x,lambda,eta){
  tmp = log(lambda) + (x-1)*log(lambda+eta*x) -lgamma(x+1) - lambda - eta*x
  return(exp(tmp))
}
# # GenPoisson density from VGAM package--I need to check if this is right
#
# PDFdgenpois <- function(x, lambda = 0, theta, log = FALSE) {
#   if (!is.logical(log.arg <- log) || length(log) != 1)
#     stop("bad input for argument 'log'")
#   rm(log)
#
#   LLL <- max(length(x), length(lambda), length(theta))
#   if (length(x)      != LLL) x      <- rep_len(x,      LLL)
#   if (length(lambda) != LLL) lambda <- rep_len(lambda, LLL)
#   if (length(theta)  != LLL) theta  <- rep_len(theta,  LLL)
#
#   llans <- -x*lambda - theta + (x-1) * log(theta + x*lambda) +
#     log(theta) - lgamma(x+1)
#   llans[x < 0] <- log(0)
#   llans[x != round(x)] <- log(0)  # x should be integer-valued
#   llans[lambda > 1] <- NaN
#   if (any(ind1 <- (lambda < 0))) {
#     epsilon <- 1.0e-9  # Needed to handle a "<" rather than a "<=".
#     mmm <- pmax(4, floor(theta/abs(lambda) - epsilon))
#     llans[ind1 & mmm < pmax(-1, -theta/mmm)] <- NaN
#     llans[ind1 & mmm < x] <- log(0)  # probability 0, not NaN
#   }
#   if (log.arg) {
#     llans
#   } else {
#     exp(llans)
#   }
# }

