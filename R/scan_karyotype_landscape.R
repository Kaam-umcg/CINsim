#' Automated scanning of karyotype landscape
#'
#' This wrapper function will allow for simulations across a range of pMisseg and selection strengths.
#'
#' @param experiment_name The name of the experiment and the folder name in which to save the results
#' @param iterations The number of replicates to run per condition.
#' @param karyotypes A matrix with cells in rows and chromosomes in colums. Rows and colums must be labelled appropriately.
#' @param euploid_ref The euploid chromosome copy number (used for calculating levels of aneuploidy).
#' @param g The maximum number of generations to be simulated.
#' @param pMisseg The range mis-segregation probabilities per chromosome copy.
#' @param pWGD The probability of a cell undergoing whole-genome duplication instead of a normal mitosis
#' @param pMissegG A vector of generation numbers during which pMisseg is set at the desired values. Will otherwise default to 0.
#' @param fit_misseg A logical whether pMisseg should be fitness-dependent.
#' @param pDivision The base probability of cell division (acts as a threshold for permanent cell arrest when division is fitness dependent).
#' @param fit_division A logical whether pDivision should be karyotype (i.e. fitness) dependent.
#' @param copy_num_boundaries A minimum and maximum number of copy numbers allowed per chromosome before a cell is deemed unviable.
#' @param selection_mode One of three possible modes, "cn_based", "rel_copy", "davoli", for copy number-based, relative copy numbed-based or Davoli-based selection.
#' @param selection_metric A fitness metric for the cn_based, rel_copy or davoli selection measure.
#' @param min_survival_euploid The minimum probability of survival for a euploid karyotype.
#' @param max_survival_euploid The maximum probability of survival for a euploid karyotype.
#' @param max_survival Maximum survival probability for the optimal karyotype score.
#' @param interval The interval of survival probabilities (range min to max).
#' @param probability_types The probabilities over which to iterate (either one of pDivision, pMisseg, pSurvival).
#' @param chrom_weights An optional vector of chromosome weights to affect cellular fitness.
#' @param qMods A vector of length 3 to modify the calculated values for pDivison, pMisseg and pSurvival respectively.
#' @param max_monosomy The minimum number of monosomies that will cause the cell to die.
#' @param min_euploid The minimum number of euploid chromosomes that must remain before the cell dies.
#' @param down_sample The maximum size of the simulated population before down-sampling.
#' @param down_sample_frac The fraction of cells with which to continue the simulation after down-sampling.
#' @param max_num_cells The maximum number of (theoretical) cells to simulate before the simulation is terminated.
#' @return A list of relevant simulation data. KaryoSim objects are saved to file.
#' @import tidyverse
#' @author Bjorn Bakker
#' @export

scan_karyotype_landscape <- function(experiment_name = "karyotype_landscape_scan",
                                     iterations = 6,
                                     num_free_cpu = 1,
                                     karyotypes = NULL,
                                     euploid_ref = 2,
                                     g = 50,
                                     pMisseg = 10^(seq(-6, -1, length.out = 30)),
                                     pWGD = 0,
                                     pMissegG = NULL,
                                     fit_misseg = FALSE,
                                     pDivision = 1,
                                     fit_division = FALSE,
                                     copy_num_boundaries = c(1, 8),
                                     selection_mode = NULL,
                                     selection_metric = NULL,
                                     min_survival_euploid = 0.1,
                                     max_survival_euploid = 0.9,
                                     max_survival = 1,
                                     interval = 0.1,
                                     probability_types = c("pSurvival"),
                                     chrom_weights = NULL,
                                     qMods = NULL,
                                     max_monosomy = NULL,
                                     min_euploid = NULL,
                                     down_sample = 5e+04,
                                     down_sample_frac = 0.25,
                                     max_num_cells = 1e+10) {
  # apply default Mps1-based selection
  if (is.null(selection_mode)) {
    selection_mode <- "cn_based"
  }
  if (is.null(selection_metric)) {
    selection_metric <- CINsim::Mps1_X
  }

  # set up coefficient list
  coefs <- make_cinsim_coeffcients(
    selection_metric = selection_metric,
    euploid_copy = euploid_ref,
    min_survival_euploid = min_survival_euploid,
    max_survival_euploid = max_survival_euploid,
    max_survival = max_survival,
    interval = interval,
    probability_types = probability_types
  )

  # create directory to store files
  dir.create(experiment_name)
  save(coefs, file = paste0(experiment_name, "/", experiment_name, "_coefficients.RData"))
  results_dir <- paste(experiment_name, "results/", sep = "/")
  dir.create(results_dir)

  # check input
  if ("pDivision" %in% probability_types) {
    pDivision <- 0 # this may seem counterintuitive, but if karyotype-dependent division is enabled
    # pDivision is the minimum probability cells need to have before dividing
  }

  # initializing list
  sim_res_list <- list()

  # run simulations and loop through all coefficients sets
  for (coefficient_set in names(coefs)) {
    # progress message
    message("Running coefficient set ", coefficient_set)

    # run multiple replicate simulations for a range of pMissegs,
    # using the sample-specific coefficients and selection metric
    sim_data <- parallelCinsim(
      karyotypes = karyotypes,
      g = g,
      pMisseg = rep(pMisseg, each = iterations),
      pWGD = pWGD,
      pMissegG = pMissegG,
      fit_misseg = fit_misseg,
      pDivision = pDivision,
      fit_division = fit_division,
      coef = coefs[[coefficient_set]],
      copy_num_boundaries = copy_num_boundaries,
      selection_mode = selection_mode,
      selection_metric = selection_metric,
      chrom_weights = chrom_weights,
      qMods = qMods,
      max_monosomy = max_monosomy,
      min_euploid = min_euploid,
      down_sample = down_sample,
      down_sample_frac = down_sample_frac,
      max_num_cells = max_num_cells,
      collect_fitness_score = TRUE
    )

    save(sim_data, file = paste0(results_dir, "/", coefficient_set, "_sim_data.RData"))

    # make a compiled version of the data
    compiled_data <- compile_karyoSimParallel(sim_data, for_summary = FALSE)

    # get sim info for pMissegs
    pMisseg_info <- sim_data %>%
      purrr::map("sim_info") %>%
      purrr::map(function(x) {
        tibble(p_misseg = as.numeric(x["pMisseg"]))
      }) %>%
      bind_rows(.id = "sim_label")

    # determine the change in mean population score
    karyotype_scores <- compiled_data$gen_measures %>%
      select(sim_label, g, mean_population_score) %>%
      inner_join(pMisseg_info, by = "sim_label")

    # compiled results
    sim_res <- (list(
      compiled_data = compiled_data,
      karyotype_scores = karyotype_scores
    ))

    # store to list
    sim_res_list[[coefficient_set]] <- sim_res
  }
  save(sim_res_list, file = paste0(results_dir, experiment_name, "_results.RData"))
  class(sim_res_list) <- "karyotype_landscape_scan"
  return(sim_res_list)
}

#' Summarizing karyotype landscape scanning output
#'
#' This function will generate a report of the karyotype landscape scanning.
#'
#' @param karyotype_landscape_scan An object of class karyotype_landscape_scan.
#' @return Plots and objects saved to file
#' @import tidyverse
#' @author Bjorn Bakker
#' @export


summarize_karyotype_landscape_scan <- function(karyotype_landscape_scan) {
  # check userinput
  if (class(karyotype_landscape_scan) != "karyotype_landscape_scan") {
    stop("An object of class karyotype_landscape_scan must be provided")
  }

  # get the karyotype scores, summarize over p_misseg and coefficients
  karyotype_scores <- karyotype_landscape_scan %>%
    purrr::map("karyotype_scores") %>%
    bind_rows(.id = "coefficient_set") %>%
    separate(coefficient_set, into = c("pDivision", "pMisseg", "pSurvival"), sep = "_")
  baseline_score <- karyotype_scores %>%
    filter(g == 0) %>%
    .$mean_population_score %>%
    mean(na.rm = TRUE)
  karyotype_scores %>%
    mutate(distance = abs(mean_population_score - baseline_score))
}
