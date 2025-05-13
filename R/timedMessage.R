#' Start timed message
#'
#' Sets up a timed message to monitor simulation progression.
#' @param ... Text for message
#' @author Aaron Taudt

startTimedMessage <- function(...) {
  x <- paste0(..., collapse = "")
  message(x, appendLF = FALSE)
  ptm <- proc.time()
  return(ptm)
}

#' Stop a timed message
#'
#' Sets up a timed message to monitor simulation progression.
#' @param ptm A started time message
#' @author Aaron Taudt

stopTimedMessage <- function(ptm) {
  time <- proc.time() - ptm
  message(" ", round(time[3], 2), "s")
}
