#' Determine the mode in a vector of values.
#'
#' This function returns the mode (most frequent element) of a vector. R does not have an in-built function to do this.
#'
#' @param x The vector from which to extract the mode.
#' @return The mode of a vector.
#' @author Bjorn Bakker

getMode <- function(x) {
    uniqv <- unique(x)
    uniqv[which.max(tabulate(match(x, uniqv)))]
}
