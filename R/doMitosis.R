#' Perform a round of mitosis with possible mis-segregations
#'
#' This function will perform mitosis on a given karyotype, apply chromosome mis-segregation events
#' if the criteria are met, and returns a matrix with two new karytoypes.
#'
#' @param karyoMat A karyotype matrix provided as a matrix on which to apply the mitosis.
#' @param pMisseg The probability of mis-segregation per chromosome copy number.
#' @return A 2-row matrix with karyotypes from daughter cells.
#' @import reshape2
#' @author Aaron Taudt & Bjorn Bakker
#' @export


doMitosis <- function(karyoMat, pMissegs = NULL) {
  # check user input
  if (is.null(pMissegs)) {
    message("pMisseg is missing - using default of 0.0025")
    pMissegs <- rep(0.0025, nrow(karyoMat))
  }

  # recalculate weighted probabilities
  pMissegsMatrix <- matrix(pMissegs,
    byrow = TRUE, nrow = nrow(karyoMat),
    ncol = ncol(karyoMat), dimnames = dimnames(karyoMat)
  )

  # calculate pMisseg weighted by copy number
  pChromWeighted <- 1 - (1 - pMissegsMatrix)^(2 * karyoMat)

  # determine missegregation event per chromosome
  gainLoss <- array(runif(length(pChromWeighted)) < pChromWeighted,
    dim = dim(karyoMat), dimnames = dimnames(karyoMat)
  )

  # summarize frequency of mis-segregation events
  missegFreq <- as.data.frame(table(rowSums(gainLoss)))
  colnames(missegFreq) <- c("numMisseg", "freq")
  chromFreq <- colSums(gainLoss) / nrow(gainLoss)

  # apply missegregations and setup daughter mat
  gainLossMat <- matrix(0, ncol = ncol(gainLoss), nrow = nrow(gainLoss))
  dimnames(gainLossMat) <- dimnames(gainLoss)
  gainLossMat[which(gainLoss)] <- replicate(n = sum(gainLoss), expr = sample(c(-1, 1), 1))

  daughterMat <- matrix(NA,
    nrow = nrow(karyoMat) * 2, ncol = ncol(karyoMat),
    dimnames = list(
      rep(rownames(karyoMat), 2),
      colnames(karyoMat)
    )
  )

  daughterMat[1:nrow(karyoMat), ] <- karyoMat + as.numeric(gainLossMat)
  daughterMat[(nrow(karyoMat) + 1):nrow(daughterMat), ] <- karyoMat - as.numeric(gainLossMat)

  # return ouput
  mitosisRes <- list(daughterMat, missegFreq, chromFreq)
  names(mitosisRes) <- c("daughterMat", "missegFreq", "chromFreq")
  return(mitosisRes)
}
