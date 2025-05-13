#' Append CINsim
#'
#' This functions will enable appending two CINsim objects to be joined together. This would allow
#' for instance continuing simulations with a karyotype matrix from a previous simulation.
#'
#' @param karyoSim1 A first karyoSim object.
#' @param karyoSim2 A second karyoSim object.
#' @return A joined karyoSim object
#' @author Bjorn Bakker
#' @import tidyverse
#' @export

append_cinsim <- function(karyoSim1 = NULL, karyoSim2 = NULL) {
  # check user input
  if (is.null(karyoSim1) | is.null(karyoSim2)) {
    stop("Two karyoSim objects are required")
  } else if (class(karyoSim1) != class(karyoSim2)) {
    stop("Objects have to be of the same class to be joined")
  } else if (ncol(karyoSim1$karyotypes) != ncol(karyoSim2$karyotypes)) {
    stop("Karyotype matrices have incompatible number of chromosome sets")
  }

  # whether simplified karyoSims are provided
  if (class(karyoSim1) == "karyoSim_simple") {
    simplified <- TRUE
  } else {
    simplified <- FALSE
  }

  # find final g for both objects
  final_g1 <- as.numeric(karyoSim1$sim_info["g"])
  final_g2 <- as.numeric(karyoSim2$sim_info["g"])
  new_g2 <- (final_g1 + 1):(final_g1 + final_g2)

  if (simplified) {
    gen_measures <- bind_rows(
      karyoSim1$gen_measures,
      karyoSim2$gen_measures %>%
        filter(g != 0) %>%
        mutate(g = new_g2)
    )

    # compile into a single object
    karyoSim <- list(
      karyotypes = karyoSim2$karyotypes,
      gen_measures = gen_measures,
      sim_info = karyoSim2$sim_info,
      selection_metric = karyoSim2$selection_metric,
      qMod = karyoSim2$qMod
    )
    class(karyoSim) <- "karyoSim_simple"
  } else {
    pop_measures <- bind_rows(
      karyoSim1$pop_measures,
      karyoSim2$pop_measures %>%
        filter(g != 0) %>%
        mutate(g = new_g2)
    )

    gen_measures <- bind_rows(
      karyoSim1$gen_measures,
      karyoSim2$gen_measures %>%
        filter(g != 0) %>%
        mutate(g = new_g2)
    )

    misseg_freq <- bind_rows(
      karyoSim1$misseg_freq,
      karyoSim2$misseg_freq %>%
        filter(g != 0) %>%
        mutate(g = new_g2)
    )

    # compile into a single object
    karyoSim <- list(
      karyotypes = karyoSim2$karyotypes,
      pop_measures = pop_measures,
      gen_measures = gen_measures,
      misseg_freq = misseg_freq,
      sim_info = karyoSim2$sim_info,
      selection_metric = karyoSim2$selection_metric,
      chrom_weights = karyoSim2$chrom_weights,
      qMods = karyoSim2$qMods
    )
    class(karyoSim) <- "karyoSim"
  }

  # return
  return(karyoSim)
}
