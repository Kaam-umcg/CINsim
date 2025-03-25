#' Convert a karyotype score into a probability of survival
#'
#' This functions calculate the survival probability given a karyotype score according to specified coefficients.
#'
#' @param karyotyepes A matrix of karyotypes with chromosomes in columns and cells in rows.
#' @param selection_mode One of three possible modes, "cn_based", "rel_copy", "davoli", for copy number-based, relative copy numbed-based or Davoli-based selection.
#' @param selection_metric A fitness metric for the cn_based, rel_copy or davoli selection measure.
#' @param chrom_weights An optional vector of chromosome weights to affect cellular fitness.
#' @return A fitness score.
#' @author Bjorn Bakker

get_score <- function(karyotypes = NULL,
                      selection_mode = NULL,
                      selection_metric = NULL,
                      chrom_weights = NULL) {

  if(is.null(karyotypes) | is.null(selection_mode) | is.null(selection_metric)) {
    stop("Selection mode or metric required")
  }

  ### cn_based selection
  if(selection_mode == "cn_based" & !is.null(selection_metric)) {

    # set-up a collection array
    probs <- array(NA, dim = dim(karyotypes), dimnames = dimnames(karyotypes))
    for(i in colnames(karyotypes)) {
      probs[, i] <- selection_metric[cbind(as.character(karyotypes[, i, drop = FALSE]), as.character(i))]
    }
    if(!is.null(chrom_weights)) {
      probs <- sweep(probs, 2, STATS = chrom_weights, FUN = "*")
    }
    sur_prob <- rowSums(probs)

    return(sur_prob)

  }

  ### relative copy number
  if(selection_mode == "rel_copy" & !is.null(selection_metric)) {

    # get median copy per cell and calculate genome-wide deviation from median copy
    chromosomes <- names(selection_metric)
    deviation <- array(NA, dim = dim(karyotypes[, chromosomes, drop = FALSE]),
                       dimnames = dimnames(karyotypes[, chromosomes, drop = FALSE]))
    deviation <- apply(karyotypes[, chromosomes, drop = FALSE], 1, function(x) {
      median_copy <- median(x)
      x <- x - median_copy
      x <- abs(x - (selection_metric*floor(median_copy/2)))
      if(!is.null(chrom_weights)) {
        x <- x * chrom_weights
      }
      x <- sum(x)
      return(x)
    })
    return(deviation)

  }

  ### davoli copy number
  if(selection_mode == "davoli" & !is.null(selection_metric)) {

    # check length of chrom_scores vector and trim if necessary
    if(length(selection_metric) != ncol(karyotypes)) {
      message("Chrom scores vector is longer than karyotype size - trimming")
      selection_metric <- selection_metric[1:ncol(karyotypes)]
    }

    # modify chromosome scores if weights are supplied
    if(!is.null(chrom_weights)) {
      selection_metric <- selection_metric * chrom_weights
    }

    # make a chromosome_score matrix
    chrom_score_mat <- matrix(data = selection_metric, ncol = ncol(karyotypes),
                              nrow = nrow(karyotypes),
                              byrow = TRUE, dimnames = dimnames(karyotypes))

    # calculate score
    score <- rowSums(karyotypes * chrom_score_mat)
    return(score)

  }

}


#' Convert a karyotype score into a probability of survival
#'
#' This functions calculate the survival probability given a karyotype score according to specified coefficients.
#'
#' @param score A karyotype score.
#' @param a A slope coefficient.
#' @param b An intercept coefficient.
#' @return A survival probability.
#' @author Bjorn Bakker

score2prob <- function(score, a, b, linear = TRUE) {
  if(linear) {
    p_sur <- a*score + b
  } else {
    p_sur <- exp((a * score) + b)/100
  }
  return(p_sur)
}


#' Get coefficients to determine karyotype survival
#'
#' This function will generate linear fit coefficients (slope, intercept) to convert
#' a karyotype fitness score into a survival probability (cn_based only).
#'
#' @param selection_metric The copy number frequency matrix.
#' @param euploid_copy The number of copies corresponding to the euploid state.
#' @param euploid_survival The probability of survival for a euploid karyotype.
#' @param max_survival The probability of survival for the most optimal karyotype.
#' @return A named vector with slope and intercept coefficients.
#' @import tidyverse
#' @author Bjorn Bakker

get_coefficients <- function(selection_metric = NULL, euploid_copy = 2,
                             euploid_survival = 0.9, max_survival = 1) {

  # check user input
  if(is.null(selection_metric)) {
    message("No custom selection metric provided - using Mps1-based default")
    selection_metric <- Mps1_X
  }

  # get euploid and optimal (max) fitness
  fitness <- tibble(p_sur = c(euploid_survival, max_survival),
                    score = c(selection_metric[as.character(euploid_copy), ] %>% sum(),
                              apply(X = selection_metric, MARGIN = 2, FUN = "max") %>% sum()))

  # alter the value of the maximum score if euploid karyotype yields maximum score already
  if(fitness$score[1] == fitness$score[2]) {
    fitness$score[1] <- fitness$score[2] * euploid_survival
  }

  # get linear fit and coefficients
  fit <- lm(data = fitness, formula = p_sur ~ score)
  fit_coef <- coef(fit)

  # formalize the coefficients - a = slope, b = intercept
  fit_coef <- signif(c(fit_coef["score"], fit_coef["(Intercept)"]), 4)
  names(fit_coef) <- c("a", "b")

  # return
  return(fit_coef)

}

#' Get list of coefficients to determine karyotype survival
#'
#' This function will generate linear fit coefficients (slope, intercept) to convert
#' a karyotype fitness score into a survival probability (cn_based only).
#'
#' @param selection_metric The copy number frequency matrix.
#' @param euploid_copy The number of copies corresponding to the euploid state.
#' @param min_survival_euploid The minimum probability of survival for a euploid karyotype.
#' @param max_survival_euploid The maximum probability of survival for a euploid karyotype.
#' @param max_survival The maximum probability of survival for the optimal karyotype.
#' @param interval The interval of the survival probabilities (range min to max).
#' @param probability_types The probabilities over which to iterate (either one of pDivision, pMisseg, pSurvival)
#' @return A list of lists with with slope and intercept coefficients.
#' @import tidyverse
#' @author Bjorn Bakker
#' @export

make_cinsim_coeffcients <- function(selection_metric = NULL, euploid_copy = 2,
                                    min_survival_euploid = 0.9, max_survival_euploid = 0.9,
                                    max_survival = 1, interval = 0.1,
                                    probability_types = c("pSurvival")) {

  # check user input
  if(is.null(selection_metric)) {
    message("No custom selection metric provided - using Mps1-based default")
    selection_metric <- Mps1_X
  }

  if(!any(probability_types %in% c("pDivision", "pMisseg", "pSurvival"))) {
    message("Given probability does not exist - defaulting to pSurvival only")
    probability_types <- c("pSurvival")
  }

  # p ranges
  coefs <- list()
  euploid_p <- seq(from = min_survival_euploid, to = max_survival_euploid, by = interval)
  for(i in 1:length(euploid_p)) {
    coefs[[i]] <- get_coefficients(selection_metric = selection_metric,
                                   euploid_copy = euploid_copy,
                                   euploid_survival = euploid_p[i],
                                   max_survival = max_survival)
  }
  names(coefs) <- paste("euploid_p", euploid_p, sep = "_")

  # if length of euploid_p is 1, return the default coefficients
  if(length(coefs) == 1) {
    coefs <- list(pDivision = coefs[[1]],
                  pMisseg = coefs[[1]],
                  pSurvival = coefs[[1]])
    return(coefs)
  } else {

    # collapse coefs into a single data frame
    coefs_df <- as.data.frame(coefs) %>% t()
    coef_length <- nrow(coefs_df)

    # make combinations of coefs
    if("pDivision" %in% probability_types) {
      pDivision <- seq(1, coef_length, 1)
    } else {
      pDivision <- coef_length
    }
    if("pMisseg" %in% probability_types) {
      pMisseg <- seq(1, coef_length, 1)
    } else {
      pMisseg <- coef_length
    }
    if("pSurvival" %in% probability_types) {
      pSurvival <- seq(1, coef_length, 1)
    } else {
      pSurvival <- coef_length
    }
    coef_combinations <- expand.grid(pDivision = pDivision,
                                     pMisseg = pMisseg,
                                     pSurvival = pSurvival) %>%
      as.matrix()

    # compile list of coefs
    coef_list <- list()
    for(i in 1:nrow(coef_combinations)) {
      coef_combination_name <- paste(euploid_p[coef_combinations[i, ]], collapse = "_")
      coef_sublist <- list(pDivision = coefs_df[coef_combinations[i, "pDivision"], ],
                           pMisseg = coefs_df[coef_combinations[i, "pMisseg"], ],
                           pSurvival = coefs_df[coef_combinations[i, "pSurvival"], ])
      coef_list[[coef_combination_name]] <- coef_sublist
    }

    return(coef_list)

  }


}
