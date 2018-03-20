#  Autocovariance function (ACVF) of FARIMA(0,d,0)
#'
#' @param d fractional differencing parameter
#' @param h lag to calculate ACVF at
#'
#' @return lag h ACVF of FARIMA(0,d,0)
#' @export
#'
#' @examples acf.farima0d0(1/4, 2)




acf.farima0d0 = function(d,h){
  if(h==0) return(1)
  k = 1:h
  prod((k-1+d)/(k-d))
}
