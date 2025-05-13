#' Simulate chromosomal instability (simplified).
#'
#' This the the main function for the CINsim package, that will allow for simulation of cell divisions and chromosomal instability.
#'
#' @param karyotypes A matrix with cells in rows and chromosomes in colums. Rows and colums must be labelled appropriately.
#' @param euploid_ref The euploid chromosome copy number (used for calculating levels of aneuploidy).
#' @param g The maximum number of generations to be simulated.
#' @param pMisseg The mis-segregation probability per chromosome copy.
#' @param pWGD The probability of a cell undergoing whole-genome duplication instead of a normal mitosis
#' @param copy_num_boundaries A minimum and maximum number of copy numbers allowed per chromosome before a cell is deemed unviable.
#' @param selection_mode One of three possible modes, "cn_based", "rel_copy", "davoli", for copy number-based, relative copy numbed-based or Davoli-based selection.
#' @param selection_metric A fitness metric for the cn_based, rel_copy or davoli selection measure.
#' @param coef A vector of length three containing the selection coefficients (slope, intercept and constant).
#' @param qMod A value to modify the calculated value for pSurvival.
#' @param down_sample The maximum size of the simulated population before down-sampling.
#' @param down_sample_frac The fraction of cells with which to continue the simulation after down-sampling.
#' @param max_num_cells The maximum number of (theoretical) cells to simulate before the simulation is terminated.
#' @param collect_fitness_score A logical whether to collect the fitness scores over time.
#' @return A karyoSim object containing all relevant information of the simulation.
#' @author Bjorn Bakker
#' @export

Cinsim_simple <- function(karyotypes = NULL,
                          g = 12,
                          pMisseg = 0.0025,
                          pWGD = 0,
                          copy_num_boundaries = c(1, 8),
                          selection_mode = NULL,
                          selection_metric = NULL,
                          coef = NULL,
                          qMod = NULL,
                          down_sample = 5e+04,
                          down_sample_frac = 0.25,
                          max_num_cells = 2e+09,
                          collect_fitness_score = FALSE) {
  # start timed message
  message("==> Starting simulation <==")
  time0 <- proc.time()

  # start timed message
  ptm <- startTimedMessage("Initializing simulation ...")

  # check user input for karyotypes
  if (is.null(karyotypes)) {
    # generate a mouse karyotype if none is provided
    message("No karyotype provided - using default female mouse karyotype")
    karyotypes <- makeKaryotypes(n = 1)
  }

  # initialize parameters
  down_sample_index <- 0

  # set-up collection lists for parameters
  true_cell_count_tally <- vector(length = g + 1)
  true_cell_count_tally[1] <- nrow(karyotypes)

  fraction_surviving <- vector(length = g + 1)
  fraction_surviving[1] <- 1

  mean_num_chr <- vector(length = g + 1)
  mean_num_chr[1] <- mean(rowSums(karyotypes), na.rm = TRUE)

  mean_population_score <- vector(length = g + 1)
  if (!is.null(selection_metric) & !is.null(selection_mode) & collect_fitness_score) {
    mean_population_score[1] <- mean(get_score(karyotypes,
      selection_mode = selection_mode,
      selection_metric = selection_metric
    ))
  } else {
    mean_population_score[1] <- "not calculated"
  }


  # check whether qMods are defined
  if (is.null(qMod)) {
    qMod <- 1
  }

  # check coefficients
  if (is.null(coef)) {
    a <- NULL
    b <- NULL
    c <- NULL
  }
  if (!is.null(selection_mode) & !is.null(selection_metric)) {
    coefs <- default_coefficients %>%
      filter(method == selection_mode)
    a <- coefs$a
    b <- coefs$b
    c <- coefs$c
    if (!is.null(coef)) {
      if (length(coef) == 3) {
        if (is.null(names(coef))) {
          names(coef) <- c("a", "b", "c")
        }
        a <- coef["a"]
        b <- coef["b"]
        c <- coef["c"]
      } else {
        stop("Coef should be a vector of length 3")
      }
    }
  }

  stopTimedMessage(ptm)

  # check max num cells
  if (is.null(max_num_cells)) {
    message("Caution: no maximum cell limit defined - simulation will halt at max g")
    max_num_cells <- "infinite"
  }

  # run simulation for g number of generations
  for (j in 1:g) {
    # start timed message
    ptm <- startTimedMessage("Processing generation ", j, " out of ", g, " ...")

    # set number of cells at start of generation
    num_cells <- nrow(karyotypes)

    # subset cells that are selected for WGD instead of a normal mitosis, if pWGD > 0
    if (pWGD > 0) {
      wgd_cells <- runif(n = num_cells) < pWGD
      karyotypes_WGD <- 2 * karyotypes[which(wgd_cells), , drop = FALSE]
      karyotypes <- karyotypes[which(!wgd_cells), , drop = FALSE]
    } else {
      karyotypes_WGD <- NULL
    }

    mitosis_res <- doMitosis(karyotypes, rep(pMisseg, num_cells))
    karyotypes <- rbind(mitosis_res$daughterMat, karyotypes_WGD)

    # determine number of daughters after mitosis
    num_daughters <- nrow(karyotypes)

    # identify viable and unviable cells
    viable_cells <- karyotype_survival_simple(
      karyotypes = karyotypes,
      copy_num_boundaries = copy_num_boundaries,
      selection_mode = selection_mode,
      selection_metric = selection_metric,
      a = a, b = b, c = c,
      qMod = qMod
    )

    # count the number of viable cells
    num_surviving_cells <- sum(viable_cells)

    # break loop if no surviving cells exist
    if (num_surviving_cells == 0) {
      message("No more surviving cells - exiting simulation")
      j <- j - 1
      break
    } else {
      # sub-select the viable cells
      karyotypes <- karyotypes[which(viable_cells), , drop = FALSE]

      # down sample if the cell count exceeds the predetermined threshold
      if (!is.null(down_sample)) {
        if (num_surviving_cells > down_sample) {
          num_sampled <- round(down_sample_frac * num_surviving_cells)
          sample_subset <- sample(1:num_surviving_cells, size = num_sampled)
          karyotypes <- karyotypes[sample_subset, , drop = FALSE]
          down_sample_index <- down_sample_index + 1
          true_cell_count <- round(num_sampled * (1 / down_sample_frac)^(down_sample_index), 0)
        } else {
          if (down_sample_index > 0) {
            true_cell_count <- round(num_surviving_cells * (1 / down_sample_frac)^(down_sample_index), 0)
          } else {
            true_cell_count <- num_surviving_cells
          }
        }
      } else {
        true_cell_count <- num_surviving_cells
      }

      # collect and store variables
      true_cell_count_tally[j + 1] <- true_cell_count
      fraction_surviving[j + 1] <- num_surviving_cells / num_daughters
      mean_num_chr[j + 1] <- mean(rowSums(karyotypes), na.rm = TRUE)
      if (!is.null(selection_mode) & !is.null(selection_metric) & collect_fitness_score) {
        mean_population_score[j + 1] <- mean(get_score(karyotypes,
          selection_mode = selection_mode,
          selection_metric = selection_metric
        ))
      } else {
        mean_population_score[j + 1] <- "not calculated"
      }
    }

    stopTimedMessage(ptm)

    # break if maximum number of cells is reached
    if (!is.null(max_num_cells)) {
      if (true_cell_count > max_num_cells) {
        message("True cell count greater than threshold: ", signif(true_cell_count, 3), " - exiting simulation")
        break
      }
    }
  }

  # setting up timed message for parameter compiline
  ptm <- startTimedMessage("Compiling all simulation parameters ...")

  # sub-select all relevant information from generations
  gen_measures <- tibble(
    g = 0:j,
    true_cell_count = as.numeric(true_cell_count_tally[1:(j + 1)]),
    fraction_surviving = fraction_surviving[1:(j + 1)],
    mean_num_chr = mean_num_chr[1:(j + 1)],
    mean_population_score = mean_population_score[1:(j + 1)]
  )

  # compile simulation info
  sim_info <- c(
    max_g = g,
    g = j,
    pMisseg = pMisseg,
    pWGD = pWGD,
    copy_num_boundaries_low = copy_num_boundaries[1],
    copy_num_boundaries_high = copy_num_boundaries[2],
    selection_mode = selection_mode,
    down_sample = down_sample,
    down_sample_frac = down_sample_frac,
    max_num_cells = max_num_cells
  )

  # make karyoSim object
  karyoSim <- list(
    karyotypes = karyotypes,
    gen_measures = gen_measures,
    sim_info = sim_info,
    selection_metric = selection_metric,
    qMod = qMod,
    coefficients = c(a = a, b = b, c = c)
  )

  # set object class
  class(karyoSim) <- "karyoSim_simple"

  stopTimedMessage(ptm)

  time1 <- proc.time() - time0
  message("==| Simulation complete - final time: ", round(time1[3], 2), "s |==")

  # return final output as karyoSim object
  return(karyoSim)
}

#' Run multiple simulations in parallel
#'
#' This function will enable parallel simulations to be run to quickly generate replicates for particular conditions.
#'
#' @param iterations The number of replicate iterations to generate.
#' @param cores The number of cores to use for the parallel simulations
#' @param karyotypes A single matrix or list of matrices with cells in rows and chromosomes in colums. Rows and colums must be labelled appropriately.
#' @param g The maximum number of generations to be simulated.
#' @param pMisseg The mis-segregation probability per chromosome copy.
#' @param pWGD The probability of a cell undergoing whole-genome duplication instead of a normal mitosis
#' @param copy_num_boundaries A minimum and maximum number of copy numbers allowed per chromosome before a cell is deemed unviable.
#' @param selection_mode One of three possible modes, "cn_based", "rel_copy", "davoli", for copy number-based, relative copy numbed-based or Davoli-based selection.
#' @param selection_metric A fitness metric for the cn_based, rel_copy or davoli selection measure.
#' @param coef A vector of length three containing the selection coefficients (slope and intercept and constant).
#' @param qMod A value modify the calculated value for pSurvival.
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

parallelCinsim_simple <- function(iterations = 6,
                                  cores = 2,
                                  karyotypes = NULL,
                                  g = 12,
                                  pMisseg = 0.0025,
                                  pWGD = 0,
                                  copy_num_boundaries = c(1, 8),
                                  selection_mode = NULL,
                                  selection_metric = NULL,
                                  coef = NULL,
                                  qMod = NULL,
                                  down_sample = 5e+04,
                                  down_sample_frac = 0.25,
                                  max_num_cells = 2e+09,
                                  collect_fitness_score = FALSE) {
  # start timed message
  message("==> Starting simulation(s) <==")
  time0 <- proc.time()

  # assess lengths of karyotypes, pmisseg and qMods in case of multiple values
  listable_value <- tibble(
    value = c("karyotypes", "pMisseg", "qMod"),
    multiple = c(
      is.list(karyotypes),
      length(pMisseg) > 1,
      length(qMod) > 1
    ),
    length = c(
      length(karyotypes),
      length(pMisseg),
      length(qMod)
    )
  )

  # adjust listable parameters if necessary
  if (any(listable_value$multiple)) {
    iterations_max <- listable_value %>%
      filter(multiple) %>%
      .$length %>%
      max()
    value_max <- listable_value %>%
      filter(multiple) %>%
      .$value

    # extend (if necessary) values to fit max iterations
    if (value_max == "karyotypes") {
      pMisseg <- rep(pMisseg, iterations_max)
      qMod <- rep(qMod, iterations_max)
    } else if (value_max == "pMisseg") {
      karyotypes <- rep(list(karyotypes), iterations_max)
      qMod <- rep(qMod, iterations_max)
    } else if (value_max == "qMod") {
      karyotypes <- rep(list(karyotypes), iterations_max)
      pMisseg <- rep(pMisseg, iterations_max)
    }

    # redefine iterations
    iterations <- iterations_max
  } else {
    karyotypes <- rep(list(karyotypes), iterations)
    pMisseg <- rep(pMisseg, iterations)
    qMod <- rep(qMod, iterations)
  }

  # reassigns cores argument due to some legacy code fixing
  num_cores <- cores

  # set-up multi-cores
  cl <- snow::makeCluster(num_cores)
  doSNOW::registerDoSNOW(cl)

  # progress bar
  pb <- txtProgressBar(min = 0, max = iterations, initial = 0, style = 3)
  progress <- function(n) {
    setTxtProgressBar(pb, n)
  }
  opts <- list(progress = progress)

  # generation list simulations
  simList <- foreach(
    iteration = 1:iterations,
    .export = c("Cinsim"),
    .options.snow = opts,
    .packages = "CINsim",
    # .verbose = TRUE,
    .inorder = FALSE
  ) %dopar% {
    Cinsim_simple(
      karyotypes = karyotypes[[iteration]],
      g = g,
      pMisseg = pMisseg[iteration],
      pWGD = pWGD,
      copy_num_boundaries = copy_num_boundaries,
      selection_mode = selection_mode,
      selection_metric = selection_metric,
      coef = coef,
      qMod = qMod[iteration],
      down_sample = down_sample,
      down_sample_frac = down_sample_frac,
      max_num_cells = max_num_cells,
      collect_fitness_score = collect_fitness_score
    )
  }

  # stop clusters
  close(pb)
  snow::stopCluster(cl)
  doParallel::stopImplicitCluster()

  # set simulation names
  names(simList) <- paste0("sim_", 1:iterations)

  time1 <- proc.time() - time0
  message("==| Simulation(s) complete - final time: ", round(time1[3], 2), "s |==")

  class(simList) <- "karyoSimParallel_simple"

  # return simList
  return(simList)
}


#' Determine cell survival (simplified)
#'
#' This function will determine the survival probability for a karyotype, depending on specific parameters.
#'
#' @param karyotypes A matrix with cells in rows and chromosomes in colums. Rows and colums must be labelled appropriately.
#' @param copy_num_boundaries A vector of length 2 with the minimum and maximum viable copy number state for any chromosome set.
#' @param selection_mode One of three possible modes, "cn_based", "rel_copy", "davoli", for copy number-based, relative copy numbed-based or Davoli-based selection.
#' @param selection_metric A fitness metric for the cn_based, rel_copy or davoli selection measure.
#' @param a Value of exponent for karyotype selection.
#' @param b Value of constant for karyotype selection.
#' @param c Value of second constant for karyotype selection.
#' @param qMod A modifier for the survival probability.
#' @return A vector of logicals whether the cells have survived or not.
#' @author Bjorn Bakker

karyotype_survival_simple <- function(karyotypes,
                                      copy_num_boundaries = c(1, 8),
                                      selection_mode = NULL,
                                      selection_metric = NULL,
                                      a = NULL, b = NULL, c = NULL,
                                      qMod = NULL) {
  # check user input
  if (is.null(karyotypes) | !is.matrix(karyotypes)) {
    stop("Karyotypes in matrix format should be provided")
  }

  # check for nullisomies and multisomies
  viable <- apply(karyotypes, 1, function(x) {
    !any(x < copy_num_boundaries[1] | x > copy_num_boundaries[2])
  })

  # check qMod value
  if (is.null(qMod)) {
    qMod <- 1
  }

  # apply relevant selection mode
  if (is.null(selection_mode)) {
    return(viable)
  } else {
    # sub-select viable karyotypes
    karyotypes_viable <- karyotypes[which(viable), , drop = FALSE]

    # get survival probabilities
    sur_prob <- setQmod(
      karyotype = karyotypes_viable,
      selection_mode = selection_mode,
      selection_metric = selection_metric,
      a = a, b = b, c = c,
      get_sur_probs = TRUE,
      qMod = qMod
    )

    # check survival
    survival <- sur_prob >= runif(n = nrow(karyotypes_viable), min = 0, max = 1)
  }

  # combine selection mode (if applied) and viability into a single indexed vector
  viable[which(viable)] <- survival
  return(viable)
}



#' Compile data from a list of parallel simulations (simplified)
#'
#' This functions will enable easy transformation of parallel Cinsim data into tidy data
#' that is compatible with the tidyverse function and plotting in ggplot2.
#'
#' @param karyoSim A karyoSim_simple or karyoSimParallel_simple object.
#' @param sim_label A custom label for the sets of simulations.
#' @return A list of three tibbles.
#' @author Bjorn Bakker
#' @import tidyverse
#' @export

compile_karyoSim_simple <- function(karyoSim,
                                    sim_label = NULL) {
  # check user input
  if (class(karyoSim) != "karyoSim_simple" & class(karyoSim) != "karyoSimParallel_simple") {
    stop("An object of class karyoSim_simple or karyoSimParallel_simple is required")
  }

  if (class(karyoSim) == "karyoSimParallel_simple") {
    if (is.null(sim_label)) {
      sim_label <- "sim_1"
    }
    karyoSim <- karyoSim %>%
      purrr::map("gen_measures") %>%
      purrr::map2_df(sim_label, function(x, y) {
        x %>%
          mutate(sim_label = y)
      })
  } else {
    karyoSim <- karyoSim$gen_measures
  }

  # return data
  return(karyoSim)
}

#' Plot CN frequencies (simplified)
#'
#' This functions will plot the frequencie of copy numbers.
#'
#' @param karyoSim_simple A karyoSim_simple or karyoSimParallel_simple object.
#' @return A plot of copy number frequencies.
#' @author Bjorn Bakker
#' @import tidyverse

plot_cn_simple <- function(karyoSim) {
  # check user input
  if (!class(karyoSim) %in% c("karyoSim_simple", "karyoSimParallel_simple")) {
    stop("An object of class karyoSim_simple or karyoSimParallel_simple is required.")
  }

  if (class(karyoSim) == "karyoSimParallel_simple") {
    # get karyotypes
    cn <- karyoSim %>%
      purrr::map("karyotypes") %>%
      purrr::map_df(function(x) {
        chromosomes <- colnames(x)

        karyotypes <- x %>%
          as.data.frame() %>%
          tibble::rownames_to_column("cell_id") %>%
          as_tibble()

        cn <- karyotypes %>%
          gather(chromosomes, key = "chromosome", value = "copy") %>%
          mutate(chromosome = factor(chromosome, levels = chromosomes)) %>%
          group_by(chromosome, copy) %>%
          dplyr::count() %>%
          ungroup() %>%
          group_by(chromosome) %>%
          mutate(freq = n / sum(n)) %>%
          ungroup()
      })

    # summarize the data to get mean frequencies and sds
    cn_sum <- cn %>%
      group_by(chromosome, copy) %>%
      summarize(
        mean_freq = mean(freq, na.rm = TRUE),
        sd_freq = sd(freq, na.rm = TRUE)
      ) %>%
      ungroup() %>%
      replace_na(list(
        mean_freq = 0,
        sd_freq = 0
      )) %>%
      group_by(chromosome) %>%
      mutate(
        mean_freq = mean_freq / sum(mean_freq),
        cum_mean_freq = cumsum(mean_freq),
        copy = factor(copy, levels = rev(sort(unique(cn$copy))))
      ) %>%
      ungroup()

    # plot bar graphs with sds
    p <- cn_sum %>%
      ggplot(aes(x = chromosome)) +
      geom_bar(aes(y = mean_freq, group = copy, fill = copy), stat = "identity") +
      geom_errorbar(aes(
        ymin = cum_mean_freq,
        ymax = cum_mean_freq + sd_freq,
        col = copy
      ), width = 0.5) +
      scale_y_continuous(breaks = seq(0, 1, 0.2)) +
      scale_fill_manual(values = copy_num_cols) +
      scale_colour_manual(values = copy_num_cols, guide = FALSE) +
      coord_cartesian(ylim = c(0, 1)) +
      labs(x = "Chromosome", y = "Frequency", fill = "Copy") +
      cinsim_theme()
  } else {
    # get karyotypes
    karyotypes_tibble <- karyoSim$karyotypes %>%
      as.data.frame() %>%
      tibble::rownames_to_column("cell_id") %>%
      as_tibble()

    # get copy number frequencies
    cn <- karyotypes_tibble %>%
      gather(colnames(karyoSim$karyotypes), key = "chromosome", value = "copy") %>%
      mutate(chromosome = factor(chromosome, levels = colnames(karyoSim$karyotypes))) %>%
      group_by(chromosome, copy) %>%
      dplyr::count() %>%
      ungroup() %>%
      group_by(chromosome) %>%
      mutate(
        freq = n / sum(n),
        copy = factor(copy, levels = rev(sort(unique(.$copy))))
      ) %>%
      ungroup()

    # plot
    p <- cn %>%
      ggplot(aes(x = chromosome)) +
      geom_bar(aes(y = freq, group = copy, fill = copy), stat = "identity") +
      scale_y_continuous(breaks = seq(0, 1, 0.2)) +
      scale_fill_manual(values = copy_num_cols) +
      scale_colour_manual(values = copy_num_cols, guide = FALSE) +
      coord_cartesian(ylim = c(0, 1)) +
      labs(x = "Chromosome", y = "Frequency", fill = "Copy") +
      cinsim_theme()
  }

  return(p)
}
