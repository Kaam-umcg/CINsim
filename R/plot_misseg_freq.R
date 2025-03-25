#' Plot mis-segregation frequencies (tidy)
#'
#' This function will plot the mis-segregation frequencies over time.
#'
#' @param karyoSim A karyoSim object that contains information of mis-segregation frequencies.
#' @param plot A logical wether to plot the results or return a list of data frame.
#' @return A plot of the mis-segregation frequencies
#' @import tidyverse
#' @import ggplot2
#' @import patchwork
#' @author Bjorn Bakker
#' @export

plot_misseg_freq <- function(karyoSim, plot = TRUE) {

  # get mis-segregation frequencies and modify the table
  suppressWarnings(misseg_freq <- karyoSim$misseg_freq %>%
    unnest(misseg_freq) %>%
    mutate(num_misseg = num_misseg %>% as.numeric() %>% as.factor()))

  suppressWarnings(chrom_freq <- karyoSim$misseg_freq %>%
    unnest(chrom_freq) %>%
    filter(freq > 0) %>%
    group_by(chromosome) %>%
    summarize(g = min(g)) %>%
    arrange(g) %>%
    mutate(chromosome = factor(chromosome, levels = rev(.$chromosome))))

  if(plot) {
    # plot the results
    p1 <- misseg_freq %>%
      ggplot(aes(x = g, y = freq)) +
      geom_bar(stat = "identity", aes(fill = num_misseg)) +
      scale_fill_viridis_d(direction = -1) +
      scale_x_continuous(breaks = c(1, seq(5, to = max(misseg_freq$g), by = 5))) +
      expand_limits(x = 1) +
      labs(x = "Generation", y = "Fraction", fill = "Mis-segregations",
           title = "Number of mis-segregations per mitosis") +
      cinsim_theme() +
      theme(legend.position = "bottom")

    p2 <- chrom_freq %>%
      ggplot(aes(x = chromosome, y = g)) +
      geom_segment(aes(y = 0, yend = g, xend = chromosome), col = "grey") +
      geom_hline(yintercept = 0, col = "black") +
      geom_point(aes(fill = chromosome), shape = 21, col = "black", size = 5) +
      scale_y_continuous(breaks = seq(0, karyoSim$sim_info["g"], by = 1)) +
      scale_fill_viridis_d(guide = FALSE) +
      coord_flip() +
      labs(x = "Chromosome", y = "Generation",
           title = "Generation of first mis-segregation") +
      cinsim_theme()

    p <- p1 + p2 + patchwork::plot_layout(nrow = 1)
    return(p)

  } else {

    return(list(misseg_freq = misseg_freq,
                chrom_freq = chrom_freq))

  }

}
