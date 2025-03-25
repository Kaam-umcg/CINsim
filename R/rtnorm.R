#' Restricted 'normal' distribution
#'
#' This function will generate integer values based on a restricted 'normal' distribution with a lower and upper boundary (1 to 8)
#'
#' @param n The number of values to generate.
#' @param mean The mean (defaults to 2).
#' @param sd The standard deviation (defaults to 1)
#' @param a The lower limit (default is 1)
#' @param b The upper limit (default is 8)

rtnorm <- function(n, mean = 2, sd = 1, a = 1, b = 8){
  qnorm(runif(n, pnorm(a, mean, sd), pnorm(b, mean, sd)), mean, sd)
}
