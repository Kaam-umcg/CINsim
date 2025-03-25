#' Get coefficients to determine karyotype survival
#'
#' This function will generate linear fit coefficients (slope, intercept, and constant) to convert
#' a karyotype fitness score into a survival probability
#'
#' @param selection_mode A selection method, either "cn_based" or "rel_copy".
#' @param euploid_survival The probability of survival for a euploid karyotype.
#' @param max_survival The probability of survival for the most optimal karyotype.
#' @return A named vector with slope, intercept and constant coefficients.
#' @import tidyverse
#' @author Bjorn Bakker
#' @export

get_coefficients <- function(selection_mode = NULL, euploid_survival = 0.9, max_survival = 1) {

  # check user input
  if(is.null(selection_mode)) {
    selection_mode <- "cn_based"
    message("Applying default selection mode: cn_based")
  } else if(!method %in% c("cn_based", "rel_copy")) {
    selection_mode <- "cn_based"
    message("Seleciton mode must be either cn_based or rel_copy - applying default selection mode: cn_based")
  }

  # sub-select fitness scores
  fitness <- optimal_fitness %>%
    filter(method %in% selection_mode) %>%
    mutate(p_sur = case_when(karyotype == "Euploid" ~ euploid_survival,
                             karyotype == "Optimal" ~ max_survival))

  # get linear fit and coefficients
  fit <- lm(data = fitness, formula = p_sur ~ score)
  fit_coef <- coef(fit)

  # formalize the coefficients - a = slope, b = intercept, c = constant (1)
  fit_coef <- signif(c(fit_coef["score"], fit_coef["(Intercept)"], 1), 4)
  names(fit_coef) <- c("a", "b", "c")

  # return
  return(fit_coef)

}
