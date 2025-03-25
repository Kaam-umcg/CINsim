#' Determine cell survival
#'
#' This function will determine the survival probability for a karyotype, depending on specific parameters.
#'
#' @param karyotypes A matrix with cells in rows and chromosomes in colums. Rows and colums must be labelled appropriately.
#' @param copy_num_boundaries A vector of length 2 with the minimum and maximum viable copy number state for any chromosome set.
#' @param selection_mode One of three possible modes, "cn_based", "rel_copy", "davoli", for copy number-based, relative copy numbed-based or Davoli-based selection.
#' @param selection_metric A fitness metric for the cn_based, rel_copy or davoli selection measure.
#' @param a Value of slope for karyotype selection.
#' @param b Value of intercept for karyotype selection.
#' @param chrom_weights A vector with chromosome weights to scale importance of one or more chromosome sets.
#' @param max_monosomy The maximum number of allowed chromosome sets to be monosomic.
#' @param min_euploid The minimum number of required chromosome sets to be euploid.
#' @param euploid_ref The euploid copy number state.
#' @param qMod A modifier for the survival probability.
#' @return A vector of logicals whether the cells have survived or not.
#' @export
#' @author Bjorn Bakker

karyotype_survival <- function(karyotypes,
                               copy_num_boundaries = c(1, 8),
                               selection_mode = NULL,
                               selection_metric = NULL,
                               a = NULL, b = NULL,
                               chrom_weights = NULL,
                               max_monosomy = NULL,
                               min_euploid = NULL,
                               euploid_ref = 2,
                               qMod = 1) {

  # check user input
  if(is.null(karyotypes) | !is.matrix(karyotypes)) {
    stop("Karyotypes in matrix format should be provided")
  }

  # check for nullisomies and multisomies
  within_boundary <- apply(karyotypes, 1, function(x) {
    !any(x < copy_num_boundaries[1] | x > copy_num_boundaries[2])
  })

  # check the number of monosmies
  if(is.null(max_monosomy)) {
    max_mono <- rep(TRUE, nrow(karyotypes))
  } else {
    max_mono <- apply(karyotypes, 1, function(x) {
      sum(x == 1) < max_monosomy
    })
  }

  # check the number of euploid chromosomes
  if(is.null(min_euploid)) {
    min_eupl <- rep(TRUE, nrow(karyotypes))
  } else {
    min_eupl <- apply(karyotypes, 1, function(x) {
      sum(x != euploid_ref) >= min_euploid
    })
  }

  # determine viable cells for downstream processing
  viable <- within_boundary & max_mono & min_eupl

  # apply relevant selection mode
  if(is.null(selection_mode)) {

    return(viable)

  } else {

    # sub-select viable karyotypes
    karyotypes_viable <- karyotypes[which(viable), , drop = FALSE]

    # get survival probabilities
    sur_prob <- setQmod(karyotype = karyotypes_viable,
                        chrom_weights = chrom_weights,
                        selection_mode = selection_mode,
                        selection_metric = selection_metric,
                        a = a, b = b,
                        get_sur_probs = TRUE,
                        qMod = qMod)

    # check survival
    survival <- sur_prob >= runif(n = nrow(karyotypes_viable), min = 0, max = 1)

  }

  # combine selection mode (if applied) and viability into a single indexed vector
  viable[which(viable)] <- survival
  return(viable)

}
