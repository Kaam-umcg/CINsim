#' Quantify levels of aneuploidy
#'
#' This function will quantify the degree of aneuploidy in a sample (calculated as the absolute deviation from euploid). Higher scores indicate more aneuploidy.
#'
#' @param karyoMat A karyotype matrix with cells in rows and chromosomes in columns.
#' @param numberOfCells The number of cells in the karyotype matrix.
#' @param euploidRef The reference euploid copy number. Defaults to 2 if none is provided.
#' @return A vector of aneuploidy scores per chromosome.
#' @author Bjorn Bakker

qAneuploidy <- function(karyoMat, numberOfCells = NULL, euploidRef = NULL) {
  # check user input and create alternate variables if necessary
  if (is.null(numberOfCells)) {
    numberOfCells <- nrow(karyoMat)
  }
  if (is.null(euploidRef)) {
    euploidRef <- 2
  }

  # quantify and return level of aneuploidy
  if (numberOfCells > 1) {
    aneuploidyScore <- colMeans(abs(karyoMat - euploidRef))
  } else {
    aneuploidyScore <- abs(karyoMat - euploidRef)
  }
  return(aneuploidyScore)
}

#' Quantify levels of heterogeneity
#'
#' This function will caculate the degree of heterogeneity in a sample. Higher scores indicate greater heterogeneity.
#'
#' @param karyoMat A karyotype matrix with cells in rows and chromosomes in columns.
#' @return A vector of heterogeneity scores per chromosome.
#' @author Bjorn Bakker

qHeterogeneity <- function(karyoMat) {
  # set number of cells
  numberOfCells <- nrow(karyoMat)

  # calculate absolute frequency of copy number states per chromosome
  if (numberOfCells > 1) {
    tabs <- apply(karyoMat, 2, function(x) {
      sort(table(x), decreasing = TRUE)
    })
  } else {
    tabs <- sapply(karyoMat, function(x) {
      sort(table(x), decreasing = TRUE)
    })
  }

  # calculate and return heterogeneity score
  if (is.list(tabs)) {
    heterogeneityScore <- unlist(lapply(tabs, function(x) {
      sum(x * 0:(length(x) - 1))
    })) / numberOfCells
  } else if (is.null(dim(tabs))) {
    heterogeneityScore <- sapply(tabs, function(x) {
      sum(x * 0:(length(x) - 1))
    }) / numberOfCells
  } else {
    heterogeneityScore <- apply(tabs, 2, function(x) {
      sum(x * 0:(length(x) - 1))
    }) / numberOfCells
  }
  return(heterogeneityScore)
}

#' Calculate the fraction of cells that deviate from the modal copy number
#'
#' This function will determine the fraction of cells that deviate from the model copy number of that chromosome.
#'
#' @param karyotypes A karyotype matrix with cells in rows and chromosomes in columns.
#' @return A vector of deviation from modal scores per chromosome.
#' @author Bjorn Bakker

dev_from_mode <- function(karyotypes) {
  # get modal copy per chromosome
  modal_copy <- apply(karyotypes, 2, function(x) {
    getMode(x)
  })

  # compare modal copy number with given copy number
  deviation <- sweep(karyotypes, 2, modal_copy, "!=")

  # get mean deviation
  deviation <- apply(deviation, 2, mean)

  return(deviation)
}
