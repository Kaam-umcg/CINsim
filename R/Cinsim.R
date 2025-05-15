#' Simulate chromosomal instability.
#'
#' This the the main function for the CINsim package, that will allow for simulation of cell divisions and chromosomal instability.
#'
#' @param karyotypes A matrix with cells in rows and chromosomes in colums. Rows and colums must be labelled appropriately.
#' @param euploid_ref The euploid chromosome copy number (used for calculating levels of aneuploidy).
#' @param g The maximum number of generations to be simulated.
#' @param pMisseg The mis-segregation probability per chromosome copy.
#' @param pWGD The probability of a cell undergoing whole-genome duplication instead of a normal mitosis
#' @param pMissegG A vector of generation numbers during which pMisseg is set at the desired values. Will otherwise default to 0.
#' @param fit_misseg A logical whether pMisseg should be fitness-dependent.
#' @param pDivision The base probability of cell division (acts as a threshold for permanent cell arrest when division is fitness dependent).
#' @param fit_division A logical whether pDivision should be karyotype (i.e. fitness) dependent.
#' @param copy_num_boundaries A minimum and maximum number of copy numbers allowed per chromosome before a cell is deemed unviable.
#' @param selection_mode One of three possible modes, "cn_based", "rel_copy", "davoli", for copy number-based, relative copy numbed-based or Davoli-based selection.
#' @param selection_metric A fitness metric for the cn_based, rel_copy or davoli selection measure.
#' @param coef A list of length three containing the karyotype fitness coefficients (slope and intercept).
#' @param chrom_weights An optional vector of chromosome weights to affect cellular fitness.
#' @param qMods A vector of length 3 to modify the calculated values for pDivison, pMisseg and pSurvival respectively.
#' @param max_monosomy The minimum number of monosomies that will cause the cell to die.
#' @param min_euploid The minimum number of euploid chromosomes that must remain before the cell dies.
#' @param down_sample The maximum size of the simulated population before down-sampling.
#' @param down_sample_frac The fraction of cells with which to continue the simulation after down-sampling.
#' @param max_num_cells The maximum number of (theoretical) cells to simulate before the simulation is terminated.
#' @param collect_fitness_score A logical whether to collect the fitness scores over time.
#' @param CnFS Boolean to indicate whether the CnFS should be calculated, defaults to TRUE
#' @param KMS Boolean to indicate whether the KMS should be calculated, defaults to FALSE
#' @param final_aneu_het_scores The final population measures (heterogeneity & aneuploidy) required for calculating the KMS.
#' @param verbose Degree of verbosity of CINsim, 0 means no verbosity
#' @return A karyoSim object containing all relevant information of the simulation.
#' @author Bjorn Bakker
#' @export

Cinsim <- function(karyotypes = NULL,
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
                   collect_fitness_score = FALSE,
                   CnFS = TRUE,
                   KMS = FALSE,
                   final_aneu_het_scores = NULL,
                   verbose = 1) {
  if (verbose >= 1) {
    # start timed message
    message("==> Starting simulation <==")
    time0 <- proc.time()

    # start timed message
    ptm <- startTimedMessage("Initializing simulation ...")
  }

  # check user input for pMisseg
  if (pMisseg >= 1 | pMisseg <= 0){
    message("pMisseg was an impossible value - setting to default (0.0025)")
    pMisseg <- 0.0025
  }

  # check user input for karyotypes
  if (is.null(karyotypes)) {
    # generate a mouse karyotype if none is provided
    message("No karyotype provided - using default female mouse karyotype")
    karyotypes <- makeKaryotypes(n = 1)
  }

  # initialize parameters
  if (is.null(chrom_weights)) {
    chromWeights <- rep(1, times = ncol(karyotypes))
  }
  down_sample_index <- 0

  # set-up collection lists for parameters
  pop_measures <- vector(mode = "list", length = g + 1)
  pop_measures[[1]] <- population_measures(
    karyotypes = karyotypes,
    euploid_ref = euploid_ref,
    g = 0
  )

  gen_measures <- vector(mode = "list", length = g + 1)
  gen_measures[[1]] <- tibble(
    g = 0,
    true_cell_count = nrow(karyotypes),
    fraction_surviving = 1,
    num_survivors = nrow(karyotypes),
    num_daughters = 0,
    mean_num_chr = mean(rowSums(karyotypes))
  )
  if (CnFS) {
    gen_measures[[1]] <- gen_measures[[1]] %>%
      mutate(CnFS = calc_CnFS(
        karyotypes = karyotypes,
        selection_metric = selection_metric
      ))
  } else {
    gen_measures[[1]] <- gen_measures[[1]] %>%
      mutate(CnFS = "not calculated")
  }
  if (KMS & is.null(final_aneu_het_scores)) {
    message("User wants KMS to be calculated but did not provide aneuploidy &
            heterogeneity scores required for this calculation - skipping KMS
            calculation.")
    KMS <- FALSE
  }
  if (KMS) {
    KMS_val <- calc_KMS(
      karyotypes = karyotypes,
      pop_measures = final_aneu_het_scores
    )
    gen_measures[[1]] <- gen_measures[[1]] %>%
      mutate(KMS = KMS_val)
  } else {
    gen_measures[[1]] <- gen_measures[[1]] %>%
      mutate(KMS = "not calculated")
  }

  if (!is.null(selection_mode) & !is.null(selection_metric) & collect_fitness_score) {
    gen_measures[[1]] <- gen_measures[[1]] %>%
      mutate(mean_population_score = mean(get_score(karyotypes,
        selection_mode = selection_mode,
        selection_metric = selection_metric,
        chrom_weights = chrom_weights
      )))
  } else {
    gen_measures[[1]] <- gen_measures[[1]] %>%
      mutate(mean_population_score = "not calculated")
  }

  misseg_freq <- vector(mode = "list", length = g)

  # check status of pdivision, and whether to apply asynchronous cell divisions
  # if (pDivision == 1) {
  #   if (fit_division == FALSE) {
  #     asynch <- FALSE
  #   } else {
  #     asynch <- TRUE
  #   }
  # } else {
  #   if (pDivision == 0 & fit_division == FALSE) {
  #     asynch <- TRUE
  #     pDivision <- 0.9
  #     message("pDivision must be greater than 0, applying default of 0.9")
  #   } else {
  #     asynch <- TRUE
  #   }
  # }

  # check status of pdivision, and whether to apply asynchronous cell divisions
  if (fit_division){
    asynch = TRUE
    if (pDivision != 0){
      message("fit_division was set to TRUE - ignoring pDivision")
      pDivision <- 0
    }
  }else{
    asynch = FALSE
  }

  # This feels very prone to bugs, as the conditionals
  # don't cover every case perfectly.
  # check whether qMods are defined
  if (!is.null(qMods)) {
    if (is.null(names(qMods))) {
      names(qMods) <- c("pDivision", "pMisseg", "pSurvival")
    }
  } else {
    qMods <- c(1, 1, 1)
    names(qMods) <- c("pDivision", "pMisseg", "pSurvival")
  }
  # check user coefficients
  if (!is.null(coef)) {
    if (!is.list(coef) | length(coef) != 3) {
      coef <- NULL
      message("Coefs must be a list of length three - using default coefficients")
    }
  }
  # assign default coefficients if custom coefs are incompatible
  if (is.null(coef)) {
    if (!is.null(selection_mode)) {
      coefs <- default_coefficients %>%
        filter(method == selection_mode) %>%
        select(a, b)
      coef <- list(
        pDivision = c(coefs$a, coefs$b),
        pMisseg = c(coefs$a, coefs$b),
        pSurvival = c(coefs$a, coefs$b)
      )
    } else {
      coef <- list(
        pDivision = c(a = 1, b = 1),
        pMisseg = c(a = 1, b = 1),
        pSurvival = c(a = 1, b = 1)
      )
    }
  }

  # check max num cells
  if (is.null(max_num_cells)) {
    message("Caution: no maximum cell limit defined - simulation will halt at max g")
    max_num_cells <- "infinite"
  }

  if (verbose >= 1) {
    stopTimedMessage(ptm)
  }

  # run simulation for g number of generations
  for (j in 1:g) {
    # start timed message
    if (verbose >= 1) {
      ptm <- startTimedMessage("Processing generation ", j, " out of ", g, " ...")
    }

    # check whether pMisseg is still applicable
    if (!is.null(pMissegG)) {
      if (j %in% pMissegG) {
        pMisseg <- pMisseg
      } else {
        pMisseg <- 0
      }
    }

    # set number of cells at start of generation
    num_cells <- nrow(karyotypes)

    # asynchronous cell division pipeline
    if (asynch) {
      # set pDivision as a function of karyotypic fitness
      if (fit_division) {
        divider_pdiv <- setQmod(
          karyotype = karyotypes,
          selection_mode = selection_mode,
          selection_metric = selection_metric,
          get_sur_probs = TRUE,
          a = coef$pDivision["a"], b = coef$pDivision["b"],
          qMod = qMods[["pDivision"]]
        )
        # set pdivision to 0 if it is a negative number
        divider_pdiv <- ifelse(divider_pdiv <= 0, 0, divider_pdiv)
        # set pdivision to 0 if pdiv is below threshold
        # the following line makes no sense for pDivision = 1, which
        # is the standard case. Removing it
        # divider_pdiv <- ifelse(divider_pdiv >= pDivision, divider_pdiv, 0)
        dividers <- runif(n = num_cells) < divider_pdiv
      } else {
        dividers <- runif(n = num_cells) < pDivision
      }
      karyotypes_dividers <- karyotypes[which(dividers), , drop = FALSE]
      karyotypes_non_dividers <- karyotypes[which(!dividers), , drop = FALSE]

      # apply mitoses if more than 1 dividing cell is present
      num_dividers <- sum(dividers)
      if (num_dividers > 0) {
        # subset cells that are selected for WGD instead of a normal mitosis, if pWGD > 0
        if (pWGD > 0) {
          wgd_cells <- runif(n = num_dividers) < pWGD
          karyotypes_WGD <- 2 * karyotypes_dividers[wgd_cells, , drop = FALSE]
          karyotypes_dividers <- karyotypes_dividers[!(wgd_cells), , drop = FALSE]
        } else {
          karyotypes_WGD <- NULL
        }

        # determine pMisseg when karyotype dependent
        if (fit_misseg) {
          pMissegs <- setQmod(
            karyotype = karyotypes_dividers,
            selection_mode = selection_mode,
            selection_metric = selection_metric,
            get_sur_probs = TRUE,
            a = coef$pMisseg["a"], b = coef$pMisseg["a"],
            qMod = qMods[["pMisseg"]]
          )
        } else {
          pMissegs <- rep(pMisseg, times = num_dividers)
        }

        # perform mitosis
        mitosis_res <- doMitosis(karyotypes_dividers, pMissegs)
        karyotypes <- rbind(karyotypes_non_dividers, mitosis_res$daughterMat, karyotypes_WGD)
      } else {
        # get mitosis_res if no cells divided
        mitosis_res <- list(
          daughterMat = karyotypes,
          missegFreq = data.frame(
            numMisseg = 0,
            freq = nrow(karyotypes)
          ),
          chromFreq = rep(0, times = ncol(karyotypes))
        )
        names(mitosis_res$chromFreq) <- colnames(karyotypes)
      }
    } else {
      # subset cells that are selected for WGD instead of a normal mitosis, if pWGD > 0
      if (pWGD > 0) {
        wgd_cells <- runif(n = num_cells) < pWGD
        karyotypes_WGD <- 2 * karyotypes[which(wgd_cells), , drop = FALSE]
        karyotypes <- karyotypes[which(!wgd_cells), , drop = FALSE]
      } else {
        karyotypes_WGD <- NULL
      }

      # determine pMisseg when karyotype dependent
      if (fit_misseg) {
        pMissegs <- setQmod(
          karyotype = karyotypes,
          selection_mode = selection_mode,
          selection_metric = selection_metric,
          get_sur_probs = TRUE,
          a = coef$pMisseg["a"], b = coef$pMisseg["b"],
          qMod = qMods[["pMisseg"]]
        )
      } else {
        pMissegs <- rep(pMisseg, times = num_cells)
      }

      mitosis_res <- doMitosis(karyotypes, pMissegs)
      karyotypes <- rbind(mitosis_res$daughterMat, karyotypes_WGD)
    }

    # determine number of daughters after mitosis
    num_daughters <- nrow(karyotypes)

    # identify viable and unviable cells
    viable_cells <- karyotype_survival(
      karyotypes = karyotypes,
      copy_num_boundaries = copy_num_boundaries,
      selection_mode = selection_mode,
      selection_metric = selection_metric,
      a = coef$pSurvival["a"], b = coef$pSurvival["b"],
      chrom_weights = chrom_weights,
      max_monosomy = max_monosomy,
      min_euploid = min_euploid,
      euploid_ref = euploid_ref,
      qMod = qMods["pSurvival"]
    )
    # count the number of viable cells
    num_surviving_cells <- sum(viable_cells)

    # break loop if no surviving cells exist
    if (num_surviving_cells == 0) { # this line is now bugged - what changed?

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

      # calculates the CnFS & KMS if flagged
      if (CnFS) {
        if (is.null(selection_metric)) {
          message("Cannot calculate CnFS, as there is no selection metric -
                  ensure your simulation has selection to calculate the CnFS")
          break
        }

        CnFS_val <- calc_CnFS(
          karyotypes = karyotypes,
          selection_metric = selection_metric
        )
      } else {
        CnFS_val <- "not calculated"
      }

      if (KMS) {
        if (is.null(final_aneu_het_scores)) {
          message("A final_aneu_het_score wasn't provided, and the KMS can thus
                  not be calculated - exiting the simulation")
          break
        }
        KMS_val <- calc_KMS(
          karyotypes = karyotypes,
          pop_measures = final_aneu_het_scores
        )
      } else {
        KMS_val <- "not calculated"
      }

      # collect and store variables
      pop_measures[[j + 1]] <- population_measures(
        karyotypes = karyotypes,
        euploid_ref = euploid_ref,
        g = j
      )

      gen_measures[[j + 1]] <- tibble(
        g = j,
        true_cell_count = true_cell_count,
        fraction_surviving = num_surviving_cells / num_daughters,
        num_survivors = num_surviving_cells,
        num_daughters = num_daughters,
        mean_num_chr = mean(rowSums(karyotypes)),
        CnFS = CnFS_val,
        KMS = KMS_val
      )
      if (!is.null(selection_mode) & !is.null(selection_metric) & collect_fitness_score) {
        gen_measures[[j + 1]] <- gen_measures[[j + 1]] %>%
          mutate(mean_population_score = mean(get_score(karyotypes,
            selection_mode = selection_mode,
            selection_metric = selection_metric,
            chrom_weights = chrom_weights
          )))
      } else {
        gen_measures[[j + 1]] <- gen_measures[[j + 1]] %>%
          mutate(mean_population_score = "not calculated")
      }

      misseg_freq[[j]] <- tibble(
        chromosome = names(mitosis_res$chromFreq),
        freq = mitosis_res$chromFreq
      ) %>%
        mutate(chromosome = factor(chromosome, levels = .$chromosome)) %>%
        nest(chrom_freq = c(chromosome, freq)) %>%
        bind_cols(tibble(
          num_misseg = mitosis_res$missegFreq$numMisseg,
          freq = mitosis_res$missegFreq$freq
        ) %>%
          mutate(freq = freq / sum(freq)) %>%
          nest(misseg_freq = c(num_misseg, freq)) %>%
          mutate(g = j))
    }
    if (verbose >= 1) {
      stopTimedMessage(ptm)
    }

    # break if maximum number of cells is reached
    if (!is.null(max_num_cells)) {
      if (true_cell_count > max_num_cells) {
        message("True cell count greater than threshold: ", signif(true_cell_count, 3), " - exiting simulation")
        break
      }
    }
  }

  # setting up timed message for parameter compiline
  if (verbose >= 1) {
    ptm <- startTimedMessage("Compiling all simulation parameters ...")
  }

  # sub-select all relevant information from generations
  pop_measures <- bind_rows(pop_measures[1:(j + 1)])
  gen_measures <- bind_rows(gen_measures[1:(j + 1)]) %>%
    mutate(true_cell_count = as.numeric(true_cell_count))
  misseg_freq <- bind_rows(misseg_freq[1:j])

  # compile simulation info
  sim_info <- c(
    euploid_ref = euploid_ref,
    max_g = g,
    g = j,
    pMisseg = pMisseg,
    pWGD = pWGD,
    pMissegG = pMissegG,
    fit_misseg = fit_misseg,
    pDivision = pDivision,
    fit_division = fit_division,
    copy_num_boundaries_low = copy_num_boundaries[1],
    copy_num_boundaries_high = copy_num_boundaries[2],
    selection_mode = selection_mode, # if selection_mode = NULL doesn't appear in object
    max_monosomy = max_monosomy,
    min_euploid = min_euploid,
    down_sample = down_sample,
    down_sample_frac = down_sample_frac,
    max_num_cells = max_num_cells
  )

  # make karyoSim object
  karyoSim <- list(
    karyotypes = karyotypes,
    pop_measures = pop_measures,
    gen_measures = gen_measures,
    misseg_freq = misseg_freq,
    sim_info = sim_info,
    selection_metric = selection_metric,
    chrom_weights = chrom_weights,
    qMods = qMods,
    coefficients = coef
  )

  # set object class
  class(karyoSim) <- "karyoSim"

  if (verbose >= 1) {
    stopTimedMessage(ptm)

    time1 <- proc.time() - time0
    message("==| Simulation complete - final time: ", round(time1[3], 2), "s |==")
  }

  # return final output as karyoSim object
  return(karyoSim)
}
