#' Run multiple simulations in parallel
#'
#' This function will enable parallel simulations to be run to quickly generate replicates for particular conditions.
#'
#' @param iterations The number of replicate iterations to generate.
#' @param cores Amount of cores that are to be assigned to CINsim. Defaults to 2.
#' @param karyotypes A single matrix or list of matrices with cells in rows and chromosomes in colums. Rows and colums must be labelled appropriately.
#' @param euploid_ref The euploid chromosome copy number (used for calculating levels of aneuploidy).
#' @param g The maximum number of generations to be simulated.
#' @param pMisseg The mis-segregation probability per chromosome copy. Can be multiple to loop over.
#' @param pWGD The probability of a cell undergoing whole-genome duplication instead of a normal mitosis
#' @param pMissegG A vector of generation numbers during which pMisseg is set at the desired values. Will otherwise be consistent.
#' @param fit_misseg A logical whether pMisseg should be fitness-dependent.
#' @param pDivision The base probability of cell division.
#' @param fit_division A logical whether pDivision should be karyotype (i.e. fitness) dependent.
#' @param copy_num_boundaries A minimum and maximum number of copy numbers allowed per chromosome before a cell is deemed unviable.
#' @param selection_mode One of three possible modes, "cn_based", "rel_copy", "davoli", for copy number-based, relative copy numbed-based or Davoli-based selection.
#' @param selection_metric A fitness metric for the cn_based, rel_copy or davoli selection measure. Usually the observed karyotype frequencies from your data
#' @param coef coefficients for alternative values of pMisseg, pSurvival & pDivision. Can be created by the (make_cinsim_coefficients) function
#' @param chrom_weights An optional vector of chromosome weights to affect cellular fitness.
#' @param qMods A single vector or list of vectors of length 3 to modify the calculated values for pDivison, pMisseg and pSurvival respectively.
#' @param max_monosomy The minimum number of monosomies that will cause the cell to die.
#' @param min_euploid The minimum number of euploid chromosomes that must remain before the cell dies.
#' @param down_sample The maximum size of the simulated population before down-sampling.
#' @param down_sample_frac The fraction of cells with which to continue the simulation after down-sampling.
#' @param max_num_cells The maximum number of (theoretical) cells to simulate before the simulation is terminated.
#' @param collect_fitness_score A logical whether to collect the fitness scores over time.
#' @return A \code{list} with karyoSim objects.
#' @import tidyverse
#' @import parallel
#' @import doParallel
#' @import doSNOW
#' @import foreach
#' @author Bjorn Bakker, Alex van Kaam
#' @export

parallelCinsim <- function(iterations = 6,
                           cores = 2,
                           karyotypes = NULL,
                           euploid_ref = 2,
                           g = 12,
                           pMisseg = 0.0025,
                           pWGD = 0,
                           pMissegG = NULL,
                           fit_misseg = FALSE,
                           pDivision = 1,
                           fit_division = FALSE,
                           copy_num_boundaries = c(1, 8),
                           selection_mode = NULL,
                           selection_metric = NULL,
                           coef = NULL,
                           chrom_weights = NULL,
                           qMods = NULL,
                           max_monosomy = NULL,
                           min_euploid = NULL,
                           down_sample = 5e+04,
                           down_sample_frac = 0.25,
                           max_num_cells = 2e+09,
                           collect_fitness_score = FALSE) {

  # start timed message
  message("==> Starting simulation(s) <==")
  time0 <- proc.time()

  # assess lengths of karyotypes, pmisseg and qMods in case of multiple values
  listable_value <- tibble(value = c("karyotypes", "pMisseg", "qMods"),
                           multiple = c(is.list(karyotypes),
                                        length(pMisseg) > 1,
                                        is.list(qMods)),
                           length = c(length(karyotypes),
                                      length(pMisseg),
                                      length(qMods)))

  # adjust listable parameters if necessary
  if(any(listable_value$multiple)) {

    iterations_max <- listable_value %>%
      filter(multiple) %>%
      .$length %>%
      max()
    value_max <- listable_value %>%
      filter(multiple) %>%
      .$value

    # extend (if necessary) values to fit max iterations
    if(value_max == "karyotypes") {
      pMisseg <- rep(pMisseg, iterations_max)
      qMods <- rep(list(qMods), iterations_max)
    } else if(value_max == "pMisseg") {
      karyotypes <- rep(list(karyotypes), iterations_max)
      qMods <- rep(list(qMods), iterations_max)
    } else if(value_max == "qMods") {
      karyotypes <- rep(list(karyotypes), iterations_max)
      pMisseg <- rep(pMisseg, iterations_max)
    }

    # redefine iterations
    iterations <- iterations_max

  } else {
    karyotypes <- rep(list(karyotypes), iterations)
    pMisseg <- rep(pMisseg, iterations)
    qMods <- rep(list(qMods), iterations)
  }

  # set-up multi-cores
  cl <- snow::makeCluster(cores)
  doSNOW::registerDoSNOW(cl)

  # progress bar
  pb <- txtProgressBar(min = 0, max = iterations, initial = 0, style = 3)
  progress <- function(n) { setTxtProgressBar(pb, n) }
  opts <- list(progress = progress)

  # generation list simulations
  simList <- foreach(iteration = 1:iterations,
                     .export = c("Cinsim"),
                     .options.snow = opts,
                     .packages = 'CINsim',
                     #.verbose = TRUE,
                     .inorder = TRUE) %dopar% {
                       Cinsim(karyotypes = karyotypes[[iteration]],
                              euploid_ref = euploid_ref,
                              g = g,
                              pMisseg = pMisseg[iteration],
                              pWGD = pWGD,
                              pMissegG = pMissegG,
                              fit_misseg = fit_misseg,
                              pDivision = pDivision,
                              fit_division = fit_division,
                              copy_num_boundaries = copy_num_boundaries,
                              selection_mode = selection_mode,
                              selection_metric = selection_metric,
                              coef = coef,
                              chrom_weights = chrom_weights,
                              qMods = qMods[[iteration]],
                              max_monosomy = max_monosomy,
                              min_euploid = min_euploid,
                              down_sample = down_sample,
                              down_sample_frac = down_sample_frac,
                              max_num_cells = max_num_cells,
                              collect_fitness_score = collect_fitness_score)
                     }

  # stop clusters
  close(pb)
  snow::stopCluster(cl)
  doParallel::stopImplicitCluster()

  # set simulation names
  names(simList) <- paste0("sim_", 1:iterations)

  time1 <- proc.time() - time0
  message("==| Simulation(s) complete - final time: ", round(time1[3], 2), "s |==")

  class(simList) <- "karyoSimParallel"

  # return simList
  return(simList)

}
