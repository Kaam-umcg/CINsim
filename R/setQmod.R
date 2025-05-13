#' Calculate qMod to result in a desired value
#'
#' This function will take calculate a qMod value such that a desired value is obtained based
#' on the selection mode and metric for a specified karyotype (i.e. a euploid karyotype). Alternatively
#' it can be used to easily calculate survival probabilities for a collection of karyotypes.
#'
#' @param karyotype A karyotype provided as a vector (default is 'euploid female mouse')
#' @param desired_value The desired value of the parameter based on the provided karyotype.
#' @param chrom_weights A vector with chromosome weights to scale importance of one or more chromosome sets.
#' @param selection_mode One of three possible modes, "cn_based", "rel_copy", "davoli", for copy number-based, relative copy numbed-based or Davoli-based selection.
#' @param selection_metric A fitness metric for the cn_based, rel_copy or davoli selection measure.
#' @param a The slope coefficient for the selection paramater.
#' @param b The intercept coefficient for the selection paramater.
#' @param get_sur_probs A logical whether to return the desired qMod, or return a survival probability.
#' @param qMod A custom qMod to modify returned survival probablities.
#' @return A qMod value.
#' @author Bjorn Bakker
#' @export
#'

setQmod <- function(karyotype = NULL, desired_value = 0.9, chrom_weights = NULL,
                    selection_mode = NULL, selection_metric = NULL, a = NULL, b = NULL,
                    get_sur_probs = FALSE, qMod = 1) {
  # check user input
  if (is.null(selection_mode) | is.null(selection_metric)) {
    stop("No selection mode or metric provided")
  } else if (is.null(karyotype)) {
    message("No karyotype provided, using default euploid female mouse")
    karyotype <- makeKaryotypes(1)
  }

  if (selection_mode %in% c("cn_based", "rel_copy", "davoli") & !is.null(selection_metric)) {
    # get coefficients
    if (is.null(a) | is.null(b)) {
      message("Using default coefficients")
      coefs <- default_coefficients %>%
        filter(method == selection_mode)
      a <- coefs$a
      b <- coefs$b
    }

    # get scores
    scores <- get_score(
      karyotypes = karyotype,
      selection_mode = selection_mode,
      selection_metric = selection_metric,
      chrom_weights = chrom_weights
    )

    # calculate the probablity of survival
    if (selection_mode %in% c("cn_based", "rel_copy")) {
      sur_prob <- qMod * score2prob(scores, a, b, linear = TRUE)
    } else {
      # davoli based selection uses exponential equation
      sur_prob <- qMod * score2prob(scores, a, b, linear = FALSE)
    }

    # calculate qMod and return value
    if (get_sur_probs) {
      return(sur_prob)
    } else {
      qMod <- desired_value / sur_prob
      names(qMod) <- NULL
      return(qMod)
    }
  } else {
    message("Selection mode not recognized: please select cn_based, rel_copy, or davoli")
  }
}
