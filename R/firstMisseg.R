#' Plot degree of contribution to deviation from first mis-segregation
#'
#' This function will plot the average (-/+ SD) deviation from the modal copy number as a function of the
#' generation of first mis-segregation. For easy interpretation of the results, all provided karyoSim
#' objects should be replicates of each other.
#'
#' @param karyoSimParallel A list of karyoSim objects that contains a summary matrix for clonality frequency per generations.
#' @param plot A logical whether to plot the results or return a table.
#' @return A plot or a tibble.
#' @import tidyverse
#' @import ggplot2
#' @import ggbeeswarm
#' @author Bjorn Bakker
#' @export
#'
plotFirstMisseg <- function(karyoSimParallel, plot = TRUE) {
  # check user input
  if (class(karyoSimParallel) != "karyoSimParallel") {
    stop("An object of class karyoSimParallel must be provided")
  }

  # get generation of first mis-segregation per chromosome per simulation
  misseg_freq <- karyoSimParallel %>%
    purrr::map("misseg_freq") %>%
    purrr::map2_df(names(.), function(x, y) {
      x %>%
        unnest(chrom_freq) %>%
        mutate(sim = y)
    }) %>%
    filter(freq > 0) %>%
    mutate(chrom_sim = paste(chromosome, sim, sep = "_")) %>%
    select(g, chrom_sim) %>%
    group_by(chrom_sim) %>%
    filter(g == min(g)) %>%
    ungroup()

  # get modal deviation at the end of the simulation
  deviation <- karyoSimParallel %>%
    purrr::map("pop_measures") %>%
    purrr::map2_df(names(.), function(x, y) {
      x %>%
        unnest(measures) %>%
        mutate(sim = y) %>%
        filter(g == max(g)) %>%
        mutate(chrom_sim = paste(chromosome, sim, sep = "_")) %>%
        select(chrom_sim, deviation)
    })

  # combine all into one tibble
  first_misseg <- inner_join(misseg_freq, deviation, by = "chrom_sim") %>%
    mutate(deviation = 100 * deviation)

  # plot, otherwise return tibble
  if (plot) {
    p <- first_misseg %>%
      ggplot(aes(x = g, y = deviation)) +
      ggbeeswarm::geom_quasirandom(alpha = 0.5, col = "grey") +
      stat_summary(fun.y = "mean", geom = "line", col = "red", size = 1) +
      stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.5, size = 1) +
      stat_summary(fun.y = mean, geom = "point", size = 2) +
      scale_x_continuous(breaks = seq(1, max(first_misseg$g), by = 1)) +
      scale_y_continuous(
        breaks = seq(0, 100, by = 10),
        limits = c(0, max(first_misseg$deviation))
      ) +
      labs(
        x = "Generation", y = "Deviation from modal copy (%)",
        title = "Modal deviation relative to time of first mis-segregation"
      ) +
      cinsim_theme()

    return(p)
  } else {
    return(first_misseg)
  }
}
