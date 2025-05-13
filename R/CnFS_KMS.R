#' Calculate the final Copy number Frequency Score for a karyoSim object.
#'
#' @param sim_results A karyoSim object from CINsim functions
#' @return The Copy number Frequency Score (CnFS) for the final
#' generation for the given simulation
#' @author Bjorn Bakker, Alex van Kaam
#' @importFrom purrr pluck
#' @export

get_final_CnFS <- function(sim_results) {
  if (!class(sim_results) == "karyoSim") {
    message("Incorrect type for given object, needs to be karyoSim")
    return(FALSE)
  }

  # gets the observed (target) karyotype frequencies
  karyo_obs <- sim_results$selection_metric

  # checks whether this simulation had a selection metric - if not, returns
  # a message and 0
  if (is.null(karyo_obs)) {
    message("There was no selection metric found for the simulation, therefore
the CnFS cannot be calculated. Setting CnFS to 0.")
    return(0)
  }

  # inits the inverse score to 0, smaller reverse score = better
  inverse_CnFS <- 0

  # and the copy number frequencies of the final generation
  final_g <- max(sim_results$pop_measures$g)
  karyo_sim_final <- sim_results$pop_measures[sim_results$pop_measures$g == final_g, ]

  # calculates the CnFS (Copy number Frequency Score) for the final generation
  for (chrom in colnames(karyo_obs)) {
    # early exit for X, since we skip that chrom for the CnFS
    if (chrom == "X") {
      next
    }

    # cast to int as we need it for comparisons
    chrom <- as.integer(chrom)
    chrom_df <- karyo_sim_final[karyo_sim_final$chromosome == chrom, ]

    for (cn in rownames(karyo_obs)) {
      freq_obs <- karyo_obs[as.character(cn), as.character(chrom)]

      # test whether there is a value for cn in chrom_df,
      # if not we assign it as 0.
      if (!(as.integer(cn) %in% chrom_df$copy)) {
        freq_sim <- 0
      } else {
        freq_sim <- chrom_df[chrom_df$copy == as.integer(cn), ]$cn_freq[[1]][["freq"]]
      }

      # takes the square of difference between observed and real frequency
      chrom_cn_CnFS <- (freq_sim - freq_obs)^2
      # and adds that to the inverse CnFS term
      inverse_CnFS <- inverse_CnFS + chrom_cn_CnFS
    }
  }
  # calculates the final CnFS from the inverse
  CnFS <- 1 / inverse_CnFS
  return(CnFS)
}

#' Calculate the Copy number Frequency Score (CnFS) given a set of karyotypes
#' and a selection metric as used by CINsim
#'
#' @param karyotypes karyotypes as given from CINsim functions
#' @param selection_metric Selection metric as used in CINsim
#' @return The Copy number Frequency Score (CnFS)
#' @author Alex van Kaam
#' @export
calc_CnFS <- function(karyotypes, selection_metric) {
  # checks whether the dimensions of karyos and selection metric match
  if (!(length(karyotypes[1, ]) == length(selection_metric[1, ]))) {
    message("Your selection metric and karyotypes don't have the same
            amount of chromosomes! Cannot calculate CnFS, setting to 0")
    return(0)
  }

  # inits the inverse score to 0, smaller reverse score = better
  inverse_CnFS <- 0

  # gets the observed frequencies from the karyotypes
  freq_per_chrom <- apply(karyotypes, 2, function(col) {
    prop_table <- table(col, useNA = "no") / length(col)
    return(prop_table)
  })
  # if every chrom has a value for each CN we can skip the following steps
  # we checks for class as if every CN is present, we get a matrix
  # otherwise we get a list of lists as the dimensions don't match, in which
  # case we need to run through the following code to unpack
  if (!((c("matrix") %in% class(freq_per_chrom) &
    c("array") %in% class(freq_per_chrom)))) {
    # Ensure all CN values are represented
    all_cn_values <- paste0(min(karyotypes):max(karyotypes))

    # Create an empty matrix with all chromosomes and CN values
    chromosomes_list <- names(freq_per_chrom)

    if (is.null(chromosomes_list)) {
      chromosomes_list <- colnames(freq_per_chrom)
    }
    freq_matrix_full <- matrix(0,
      nrow = length(chromosomes_list), ncol = length(all_cn_values),
      dimnames = list(
        chromosomes_list,
        paste0(min(as.integer(all_cn_values)):max(as.integer(all_cn_values)))
      )
    )

    # Fill in the matrix with the actual frequencies from the list of tables
    for (i in 1:length(freq_per_chrom)) {
      chrom_name <- names(freq_per_chrom)[i]
      # same story for this conditional
      if (is.null(chrom_name)) {
        chrom_name <- colnames(freq_per_chrom)[i]
      }
      chrom_data <- freq_per_chrom[[i]]

      for (cn in names(chrom_data)) {
        cn_index <- which(all_cn_values == cn)
        freq_matrix_full[chrom_name, cn_index] <- chrom_data[cn]
      }
    }
    freq_per_chrom <- t(freq_matrix_full)
  }

  # calculates the CnFS (Copy number Frequency Score)
  for (chrom in colnames(selection_metric)) {
    # early exit for X, since we skip that chromosome for the CnFS
    if (chrom == "X") {
      next
    }
    for (cn in rownames(selection_metric)) {
      # gets the frequency of a CN for a chromosome from the simulation
      freq_sim <- freq_per_chrom[, chrom] %>%
        purrr::pluck(cn, .default = 0) # gives 0 if it doesn't exist/out of bounds
      # this value should always exist
      freq_obs <- selection_metric[cn, chrom]

      # takes the square of difference between observed and real frequency
      chrom_cn_CnFS <- (freq_sim - freq_obs)^2
      # and adds that to the inverse CnFS term
      inverse_CnFS <- inverse_CnFS + chrom_cn_CnFS
    }
  }
  # calculates the final CnFS from the inverse
  CnFS <- 1 / inverse_CnFS

  return(CnFS)
}

#' @title Calcualte KMS from karyotype & population measures
#' Calculate the Karyotype Matching Score (KMS) given a set of karyotypes
#' and observed aneuploidy/heterogeneity scores of the simulated population.
#'
#' @param karyotypes karyotypes as given from CINsim functions
#' @param pop_measures Selection metric as used in CINsim
#' @return The Copy number Frequency Score (CnFS)
#' @author Alex van Kaam
#' @export
calc_KMS <- function(karyotypes, pop_measures, euploid_ref = 2) {
  # gets the aneu and het scores for the given karyos
  aneu_sim <- as.vector(qAneuploidy(karyotypes, euploid_ref = euploid_ref))
  het_sim <- as.vector(qHeterogeneity(karyotypes))

  # reshapes the pop_measures object to separate vars
  aneu_pop <- pop_measures[pop_measures$parameter == "aneuploidy", ]$score
  het_pop <- pop_measures[pop_measures$parameter == "heterogeneity", ]$score

  # inits the KMS, smaller is better
  inverse_KMS <- 0

  for (chrom in seq_along(1:(length(aneu_sim) - 1))) { # -1 for skipping X
    # calculates the score for a single chrom
    inverse_KMS_score <- (aneu_pop[[as.integer(chrom)]] - aneu_sim[[as.integer(chrom)]])^2 +
      (het_pop[[as.integer(chrom)]] - het_sim[[as.integer(chrom)]])^2

    inverse_KMS <- inverse_KMS + inverse_KMS_score
  }
  KMS <- 1 / inverse_KMS
  return(KMS)
}
