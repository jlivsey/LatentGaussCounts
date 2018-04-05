temp_Hermite_poly <- function(k,x){
  # Polys is preloaded as data only visable to LatentGaussCounts package
  sum( as.numeric(Polys[[k]]) * (x^(0:(k-1))) )
}

#' Hermite Polynomials
#'
#' evaluate the k-th Hermite Polynomial at x
#'
#' @param k Hermite polynomial
#' @param x input
#'
#' @return value of H_k(x)
#' @export
#'

Hermite_poly = Vectorize(temp_Hermite_poly, vectorize.args = "x")
