#' Make a matrix of karyotypes
#'
#' This function will create a matrix of karyotypes that is compatible with the CINsim pipeline. Currently supported species are mouse (20 chromosomes) and human (23 chromsomes).
#'
#' @param numCell The desired number of cells in the matrix.
#' @param species The desired species, either mouse or human (20 or 23 chromosome pairs). Otherwise specify the number of chromosome pairs to use.
#' @param copies The number of copies for all chromosomes in the karyotype matrix. For random karyotypes this will correspond to the ploidy around which random karyotypes will generated.
#' @param random A logical whether to generate random karyotypes (the copies argument will provide the base ploidy (25 percent deviance around the base ploidy will be used).
#' @param varianceMultiplier The degree of variance (in standard deviations) in randomly generated karyotypes.
#' @return A matrix of karyotypes according to the input parameters.
#' @author Bjorn Bakker
#' @export

makeKaryotypes <- function(numCell = 20, species = "mouse", copies = 2,
                           random = FALSE, varianceMultiplier = 0.5) {

  # check user input to determine the number of chromosome pairs
  if(species == "mouse") {
    numChr <- 20
  } else if(species == "human") {
    numChr <- 23
  } else {
    numChr <- species
  }

  # generate random karyotypes if requested
  if(random) {

    # assign random karyotypes around the ploidy count
    karyotypeMat <- t(replicate(n = numCell,
                                expr = round(rtnorm(numChr, mean = copies, sd = varianceMultiplier), digits = 0)))

    # check for cells with negative copy numbers or greater than 9 and remove them; make additional karyotypes until
    # the requested number of karyotypes has been generated
    impossibleCells <- apply(karyotypeMat, 1, function(x) {any(x < 0 | x > 9)})
    numImpossibleCells <- sum(impossibleCells)
    while(numImpossibleCells > 0) {
      newCells <- t(replicate(n = numImpossibleCells,
                              expr = round(rtnorm(numChr, mean = copies, sd = varianceMultiplier), digits = 0)))
      karyotypeMat <- rbind(karyotypeMat[!impossibleCells, ], newCells)
      impossibleCells <- apply(karyotypeMat, 1, function(x) {any(x < 0 | x > 9)})
      numImpossibleCells <- sum(impossibleCells)
    }


  } else {
    karyotypes <- rep(copies, times = numCell*numChr)
    karyotypeMat <- matrix(karyotypes, ncol = numChr)
  }

  # compile the karyotype matrix
  colnames(karyotypeMat) <- c(1:(numChr - 1), "X")
  rownames(karyotypeMat) <- paste0("cell_", 1:numCell)
  return(karyotypeMat)

}
