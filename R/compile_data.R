#' Compile data from CINsim
#'
#' This functions will enable easy transformation of parallel Cinsim data into tidy data
#' that is compatible with the tidyverse functions and plotting in ggplot2.
#'
#' @param karyoSim A karyoSim object.
#' @param sim_label A custom label for the simulation.
#' @return A list of three tibbles.
#' @author Bjorn Bakker
#' @import tidyverse
#' @export

compile_karyoSim <- function(karyoSim,
                             sim_label = NULL) {
  # check user input
  if (class(karyoSim) != "karyoSim") {
    stop("An object of class karyoSim is required")
  }

  # check sim labels, otherwise use the names of the elements in list
  if (is.null(sim_label)) {
    sim_label <- "sim"
  }

  # gather per chromosome and genomewide information
  pop_measures_chrom <- karyoSim$pop_measures %>%
    unnest(measures) %>%
    gather(-g, -chromosome, key = "parameter", value = "value") %>%
    mutate(sim_label = sim_label)

  pop_measures_gw <- karyoSim$pop_measures %>%
    unnest(measures) %>%
    gather(-g, -chromosome, key = "parameter", value = "value") %>%
    group_by(g, parameter) %>%
    summarize(value = mean(value, na.rm = TRUE)) %>%
    ungroup() %>%
    remove_missing(na.rm = TRUE) %>%
    mutate(
      sim_label = sim_label,
      parameter = factor(parameter, levels = c(
        "aneuploidy",
        "heterogeneity",
        "deviation"
      ))
    )

  gen_measures <- karyoSim$gen_measures %>%
    gather(-g, key = "parameter", value = "value") %>%
    mutate(
      sim_label = sim_label,
      value = as.numeric(value)
    )

  # return list
  compiled_data <- list(
    pop_measures_chrom = pop_measures_chrom,
    pop_measures_gw = pop_measures_gw,
    gen_measures = gen_measures
  )
  return(compiled_data)
}

#' Compile data from a list of parallel simulations
#'
#' This functions will enable easy transformation of parallel Cinsim data into tidy data
#' that is compatible with the tidyverse functions and plotting in ggplot2.
#'
#' @param karyoSimParallel A karyoSimParallel object.
#' @param sim_label A custom label for the sets of simulations.
#' @param for_summary A logical whether the data will be summarized later.
#' @return A list of three tibbles.
#' @author Bjorn Bakker
#' @import tidyverse
#' @export

compile_karyoSimParallel <- function(karyoSimParallel,
                                     sim_label = NULL,
                                     for_summary = FALSE) {
  # check user input
  if (class(karyoSimParallel) != "karyoSimParallel") {
    stop("An object of class karyoSimParallel is required")
  }

  # check sim labels, otherwise use the names of the elements in list
  if (is.null(sim_label) & is.null(names(karyoSimParallel))) {
    sim_label <- paste("set", 1:length(karyoSimParallel), sep = "_")
  } else if (is.null(sim_label)) {
    sim_label <- names(karyoSimParallel)
  }

  # compile all metrics into a single list of tidy data
  pop_measures <- karyoSimParallel %>%
    purrr::map("pop_measures") %>%
    purrr::map2_df(sim_label, function(x, y) {
      x <- x %>%
        mutate(sim_label = y)
    })

  gen_measures <- karyoSimParallel %>%
    purrr::map("gen_measures") %>%
    purrr::map2_df(sim_label, function(x, y) {
      x <- x %>%
        mutate(sim_label = y)
    })

  compiled_data <- list(
    pop_measures = pop_measures,
    gen_measures = gen_measures
  )

  # return list
  if (for_summary) {
    # gather per chromosome and genomewide information
    pop_measures_chrom <- compiled_data$pop_measures %>%
      unnest(measures) %>%
      gather(-g, -chromosome, -sim_label, key = "parameter", value = "value") %>%
      select(g, chromosome, parameter, value, sim_label)

    pop_measures_gw <- compiled_data$pop_measures %>%
      unnest(measures) %>%
      gather(-g, -chromosome, -sim_label, key = "parameter", value = "value") %>%
      group_by(g, parameter, sim_label) %>%
      summarize(value = mean(value, na.rm = TRUE)) %>%
      ungroup() %>%
      # replace_na(list(mean_value = 0, sd_value = 0)) %>%
      remove_missing(na.rm = TRUE) %>%
      mutate(parameter = factor(parameter, levels = c(
        "aneuploidy",
        "heterogeneity",
        "deviation"
      ))) %>%
      select(g, parameter, value, sim_label)

    gen_measures <- compiled_data$gen_measures %>%
      gather(-g, -sim_label, key = "parameter", value = "value") %>%
      select(g, parameter, value, sim_label) %>%
      mutate(value = as.numeric(value))

    # return compiled list
    compiled_data <- list(
      pop_measures_chrom = pop_measures_chrom,
      pop_measures_gw = pop_measures_gw,
      gen_measures = gen_measures
    )
    return(compiled_data)
  } else {
    return(compiled_data)
  }
}

#' Summarize population metrics of a karyoSimParallel object
#'
#' This functions will enable easy transformation of parallel Cinsim data into tidy summarized data
#' that is compatible with the tidyverse functions and plotting in ggplot2.
#'
#' @param karyoSimParallel A karyoSimParallel objects.
#' @param sim_label A custom label for the sets of simulations.
#' @return A list of three tibbles.
#' @author Bjorn Bakker
#' @import tidyverse
#' @export

summarize_karyoSimParallel <- function(karyoSimParallel, sim_label = NULL) {
  # check user input
  if (class(karyoSimParallel) != "karyoSimParallel") {
    stop("An object of class karyoSimParallel is required")
  }

  # check sim labels, otherwise use the names of the elements in list
  if (is.null(sim_label)) {
    sim_label <- "summary"
  }

  # compile data
  compiled_data <- compile_karyoSimParallel(karyoSimParallel, for_summary = TRUE)

  # summarize generation measures
  gen_measures <- compiled_data$gen_measures %>%
    select(-sim_label) %>%
    gather(-g, key = "parameter", value = "value") %>%
    group_by(g, parameter) %>%
    summarize(
      mean_value = mean(value, na.rm = TRUE),
      sd_value = sd(value, na.rm = TRUE)
    ) %>%
    ungroup() %>%
    # replace_na(list(mean_value = 0, sd_value = 0)) %>%
    remove_missing(na.rm = TRUE) %>%
    mutate(sim_label = sim_label)

  # summarize population measures per chromosome
  pop_measures <- compiled_data$pop_measures %>%
    select(g, measures) %>%
    unnest() %>%
    gather(-g, -chromosome, key = "parameter", value = "value")

  pop_measures_chrom <- pop_measures %>%
    group_by(g, chromosome, parameter) %>%
    summarize(
      mean_value = mean(value, na.rm = TRUE),
      sd_value = sd(value, na.rm = TRUE)
    ) %>%
    ungroup() %>%
    # replace_na(list(mean_value = 0, sd_value = 0)) %>%
    remove_missing(na.rm = TRUE) %>%
    mutate(
      sim_label = sim_label,
      parameter = factor(parameter, levels = c(
        "aneuploidy",
        "heterogeneity",
        "deviation"
      ))
    )

  # summarize population measures genomewide
  pop_measures_gw <- pop_measures %>%
    group_by(g, parameter) %>%
    summarize(
      mean_value = mean(value, na.rm = TRUE),
      sd_value = sd(value, na.rm = TRUE)
    ) %>%
    ungroup() %>%
    # replace_na(list(mean_value = 0, sd_value = 0)) %>%
    remove_missing(na.rm = TRUE) %>%
    mutate(
      sim_label = sim_label,
      parameter = factor(parameter, levels = c(
        "aneuploidy",
        "heterogeneity",
        "deviation"
      ))
    )

  # compile data into a single list
  compiled_data <- list(
    gen_measures = gen_measures,
    pop_measures_chrom = pop_measures_chrom,
    pop_measures_gw = pop_measures_gw
  )

  class(compiled_data) <- "karyoSimSummary"
  return(compiled_data)
}

#' Summarize population metrics of multiple karyoSimParallel objects
#'
#' This functions will enable easy transformation of parallel Cinsim data into tidy summarized data
#' that is compatible with the tidyverse functions and plotting in ggplot2.
#'
#' @param karyoSimParallel A list of karyoSimParallel objects.
#' @param sim_label A custom label for the sets of simulations.
#' @return A list of three tibbles.
#' @author Bjorn Bakker
#' @import tidyverse
#' @export

summarize_karyoSimParallel_list <- function(karyoSimParallel, sim_label = NULL) {
  # check user input
  object_class <- purrr::map_chr(karyoSimParallel, class)
  if (any(object_class != "karyoSimParallel")) {
    stop("A list with objects of class karyoSimParallel is required")
  }

  # check sim_label names
  if (is.null(sim_label) & is.null(names(karyoSimParallel))) {
    sim_label <- paste0("sim_", 1:length(karyoSimParallel))
  } else if (is.null(sim_label)) {
    sim_label <- names(karyoSimParallel)
  } else if (length(sim_label) != length(karyoSimParallel)) {
    message("Length of sim_labels inconsistent with length of list - using default labels")
    sim_label <- paste0("sim_", 1:length(karyoSimParallel))
  }

  # summarize data and combine all data into a single big list of data frames
  summarized_data <- purrr::map2(karyoSimParallel, sim_label, function(x, y) {
    summarize_karyoSimParallel(karyoSimParallel = x, sim_label = y)
  })
  compiled_data <- list(
    gen_measures = summarized_data %>%
      purrr::map("gen_measures") %>%
      bind_rows(),
    pop_measures_chrom = summarized_data %>%
      purrr::map("pop_measures_chrom") %>%
      bind_rows(),
    pop_measures_gw = summarized_data %>%
      purrr::map("pop_measures_gw") %>%
      bind_rows()
  )

  # return compiled data
  class(compiled_data) <- "karyoSimSummary"
  return(compiled_data)
}

#' Compile the metrics of multiple karyoSimParallel objects within distinct inputs (i.e. non-replicates)
#'
#' This functions will enable easy transformation of parallel Cinsim data into tidy compiled data
#' that is compatible with the tidyverse functions and plotting in ggplot2. This function is intended for
#' parallel simulations that each have a unique parameter (e.g. they have looped over values for pMisseg).
#'
#' @param karyoSimParallel A list of karyoSimParallel objects.
#' @param sim_label A custom label for the sets of simulations.
#' @return A list of four tibbles.
#' @author Bjorn Bakker
#' @import tidyverse
#' @export

compile_iterative_karyoSimParallel <- function(karyoSimParallel, sim_label = NULL) {
  # check user input
  object_class <- purrr::map_chr(karyoSimParallel, class)
  if (any(object_class != "karyoSimParallel")) {
    stop("A list with objects of class karyoSimParallel is required")
  }

  # check sim_label names
  if (is.null(sim_label) & is.null(names(karyoSimParallel))) {
    sim_label <- paste0("rep_", 1:length(karyoSimParallel))
  } else if (is.null(sim_label)) {
    sim_label <- names(karyoSimParallel)
  } else if (length(sim_label) != length(karyoSimParallel)) {
    message("Length of sim_labels inconsistent with length of list - using default labels")
    sim_label <- paste0("rep_", 1:length(karyoSimParallel))
  }

  # collect simulation info
  sim_info <- purrr::map2_df(karyoSimParallel, sim_label, function(x1, x2) {
    x1 %>%
      purrr::map("sim_info") %>%
      purrr::map2_df(names(x1), function(y1, y2) {
        p_misseg <- y1["pMisseg"] %>% as.numeric()
        tibble(
          rep = x2,
          sim = y2,
          p_misseg = p_misseg
        )
      })
  }) %>%
    mutate(sim_label = paste(rep, sim, sep = "-"))

  # compile data
  compiled_data <- karyoSimParallel %>%
    purrr::map2(sim_label, function(x1, x2) {
      x2 <- paste(x2, names(x1), sep = "-")
      compile_karyoSimParallel(x1, sim_label = x2, for_summary = FALSE)
    }) %>%
    map(function(x) {
      x %>%
        map(inner_join, sim_info, by = "sim_label")
    })

  # cn_data
  cn_data <- karyoSimParallel %>%
    purrr::map2_df(sim_label, function(x1, x2) {
      x2 <- paste(x2, names(x1), sep = "-")
      cn_freq <- x1 %>%
        purrr::map("pop_measures") %>%
        purrr::map2_df(x2, function(y1, y2) {
          y1 %>%
            dplyr::select(g, cn_freq) %>%
            mutate(sim_label = y2)
        })
    }) %>%
    inner_join(sim_info, by = "sim_label") %>%
    unnest()

  # compile it all
  compiled_data_complete <- list(
    pop_measures_gw = compiled_data %>%
      purrr::map_df("pop_measures_gw"),
    pop_measures_chrom = compiled_data %>%
      purrr::map_df("pop_measures_chrom"),
    gen_measures = compiled_data %>%
      purrr::map_df("gen_measures"),
    cn_freq = cn_data
  ) %>%
    map(function(x) {
      x %>% select(-sim_label)
    })
  class(compiled_data_complete) <- "karyoSimParallel_iterative_summarized"
  return(compiled_data_complete)
}
