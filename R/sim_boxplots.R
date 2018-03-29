#' Simulation boxplot function
#'
#' @param data A dataframe with columns "parameter" (a character eg "lambda")
#'                                "estimation.method" (eg "Gaussian-lik")
#'                                 "value" (numeric paramter estimates)
#'
#' @return a boxplot
#' @export
#'


sim_boxplots <- function(data){
  require(ggplot2)
  ggplot(data, aes(x=parameter, y=value, fill=estimation.method)) +
    geom_boxplot()
}
