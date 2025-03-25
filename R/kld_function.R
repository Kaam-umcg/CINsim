#' Determine the karyotype landscape dissimilarity (KLD) score
#'
#' This function returns a score that quantifies the  dissimilarity of a karyotype landscape compared to observed
#' aneuploidy and heterogeneity scores (on a chromosome basis). It can be used to determine whether results of a
#' simulation yield a karyotype landscape that resembles a tumour sample.
#'
#' @param karyoSim A karyoSim object containing simulation information.
#' @param km_reference A table that contains the karyotype measures (aneuploidy and heterogeneity score) per chromosome.
#' @return A data frame with the karyotype dissimilarity index (KLD) per generation (absolute value and normalized to g0).
#' @author Bjorn Bakker

calculate_kld <- function(karyoSim, km_reference = NULL) {

  # if no karyotype measure reference is provided the default Mps1 table is used
  if(is.null(km_reference)) {
    km_reference <- CINsim::karyotype_measures
  }

  # determine the KLD (sum of absolute difference between observed and simulated karyotype measures)
  KLD <- karyoSim$pop_measures %>%
    select(g, measures) %>%
    unnest() %>%
    gather(-g, -chromosome, key = "parameter", value = "simulated_score") %>%
    filter(parameter %in% c("aneuploidy", "heterogeneity")) %>%
    inner_join(karyotype_measures, by = c("chromosome", "parameter")) %>%
    mutate(score_diff = abs(score - simulated_score)) %>%
    group_by(g) %>%
    summarize(kld = sum(score_diff))
  KLD_0 <- KLD %>%
    filter(g == 0) %>%
    .$kld
  KLD <- KLD %>%
    mutate(kld_norm = kld/KLD_0)
  return(KLD)

}
