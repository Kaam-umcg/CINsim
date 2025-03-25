#' Plot population measures
#'
#' This functions will plot the various population measures over time.
#'
#' @param karyoSim A karyoSim, karyoSimParallel, or list of karyoSimParallel objects.
#' @param sim_label A label for karyoSimParallel data.
#' @param summarize A logical whether to summarize a karyoSimParallel object - karyoSimParallel lists will always be summarized.
#' @param per_chromosome A logical whether to plot data per chromosome or genome-wide.
#' @param plot A logical whether to print the plot, or return a list of plots.
#' @return A summary statistics plot.
#' @author Bjorn Bakker
#' @import tidyverse
#' @import ggplot2
#' @import RColorBrewer
#' @export

plot_cinsim_summary <- function(karyoSim,
                                sim_label = NULL,
                                summarize = TRUE,
                                per_chromosome = TRUE,
                                plot = FALSE) {

  # check user data and compile tidysim if necessary
  if(class(karyoSim) == "karyoSim") {
    karyoSim <- compile_karyoSim(karyoSim, sim_label = sim_label)
    summarized <- FALSE
  } else if(class(karyoSim) == "karyoSimParallel") {
    if(summarize) {
      sim_label <- NULL
      karyoSim <- summarize_karyosimParallel(karyoSim, sim_label = sim_label)
      summarized <- TRUE
    } else {
      karyoSim <- compile_karyoSimParallel(karyoSim, sim_label = sim_label)
      summarized <- FALSE
    }
  } else if(is.list(karyoSim)) {
    if(all(purrr::map_chr(karyoSim, class) == "karyoSimParallel")) {
      karyoSim <- summarize_karyoSimParallel_list(karyoSim, sim_label = sim_label)
      summarized <- TRUE
    } else {
      stop("Data could not be processed correctly - data yields incorrect format")
    }
  } else {
    stop("A karyoSim, karyoSimParallel, or list of karyoSimParallel objects is required")
  }

  # factor and relevel the sim labels as given in order by the user
  if(class(karyoSim) != "karyoSim") {
    karyoSim <- karyoSim %>%
      purrr::map(function(x) {
        x %>% mutate(sim_label = factor(sim_label, levels = unique(x$sim_label)))
      })
  }

  # plot summarized plots - karyoSimParallel and lists
  if(summarized) {

    cell_count <- karyoSim$gen_measures %>%
      filter(parameter == "true_cell_count") %>%
      ggplot(aes(x = g, group = sim_label, col = sim_label)) +
      geom_ribbon(aes(ymin = mean_value - sd_value,
                      ymax = mean_value + sd_value,
                      fill = sim_label),
                  alpha = 0.25, col = NA, show.legend = FALSE) +
      geom_line(aes(y = mean_value), size = 1) +
      labs(x = "Generation", y = "Population size",
           title = "Growth", col = "Sim label") +
      scale_y_log10(breaks = scales::trans_breaks("log10", function(x) {10^x}),
                    labels = scales::trans_format("log10", scales::math_format(10^.x))) +
      cinsim_theme()

    fraction_surviving <- karyoSim$gen_measures %>%
      filter(parameter == "fraction_surviving") %>%
      ggplot(aes(x = g, group = sim_label, col = sim_label)) +
      geom_ribbon(aes(ymin = mean_value - sd_value,
                      ymax = mean_value + sd_value,
                      fill = sim_label),
                  alpha = 0.25, col = NA, show.legend = FALSE) +
      geom_line(aes(y = mean_value), size = 1) +
      labs(x = "Generation", y = "Fraction",
           title = "Survival", col = "Sim label") +
      cinsim_theme()

    mean_num_chr <- karyoSim$gen_measures %>%
      filter(parameter == "mean_num_chr") %>%
      ggplot(aes(x = g, group = sim_label, col = sim_label)) +
      geom_ribbon(aes(ymin = mean_value - sd_value,
                      ymax = mean_value + sd_value,
                      fill = sim_label),
                  alpha = 0.25, col = NA, show.legend = FALSE) +
      geom_line(aes(y = mean_value), size = 1) +
      labs(x = "Generation", y = "Number",
           title = "Chromosome count", col = "Sim label") +
      cinsim_theme()

    # plot data per chromosome
    if(per_chromosome) {

      aneuploidy <- karyoSim$pop_measures_chrom %>%
        filter(parameter == "aneuploidy") %>%
        ggplot(aes(x = g, group = chromosome, col = chromosome)) +
        geom_ribbon(aes(ymin = mean_value - sd_value,
                        ymax = mean_value + sd_value,
                        fill = chromosome),
                    alpha = 0.1, col = NA, show.legend = FALSE) +
        geom_line(aes(y = mean_value), size = 1) +
        facet_wrap(~sim_label) +
        scale_color_viridis_d() +
        scale_fill_viridis_d() +
        labs(x = "Generation", y = "Score",
             title = "Aneuploidy", col = "Chromosome") +
        cinsim_theme() +
        theme(legend.position = "bottom")

      heterogeneity <- karyoSim$pop_measures_chrom %>%
        filter(parameter == "heterogeneity") %>%
        ggplot(aes(x = g, group = chromosome, col = chromosome)) +
        geom_ribbon(aes(ymin = mean_value - sd_value,
                        ymax = mean_value + sd_value,
                        fill = chromosome),
                    alpha = 0.1, col = NA, show.legend = FALSE) +
        geom_line(aes(y = mean_value), size = 1) +
        facet_wrap(~sim_label) +
        scale_color_viridis_d() +
        scale_fill_viridis_d() +
        labs(x = "Generation", y = "Score",
             title = "Heterogeneity", col = "Chromosome") +
        cinsim_theme() +
        theme(legend.position = "bottom")

      deviation <- karyoSim$pop_measures_chrom %>%
        filter(parameter == "deviation") %>%
        ggplot(aes(x = g, group = chromosome, col = chromosome)) +
        geom_ribbon(aes(ymin = mean_value - sd_value,
                        ymax = mean_value + sd_value,
                        fill = chromosome),
                    alpha = 0.1, col = NA, show.legend = FALSE) +
        geom_line(aes(y = mean_value), size = 1) +
        facet_wrap(~sim_label) +
        scale_color_viridis_d() +
        scale_fill_viridis_d() +
        labs(x = "Generation", y = "Score",
             title = "Deviation", col = "Chromosome") +
        cinsim_theme() +
        theme(legend.position = "bottom")

      # plot data per genome-wide
    } else {

      aneuploidy <- karyoSim$pop_measures_gw %>%
        filter(parameter == "aneuploidy") %>%
        ggplot(aes(x = g, group = sim_label, col = sim_label)) +
        geom_ribbon(aes(ymin = mean_value - sd_value,
                        ymax = mean_value + sd_value,
                        fill = sim_label),
                    alpha = 0.25, col = NA, show.legend = FALSE) +
        geom_line(aes(y = mean_value), size = 1) +
        labs(x = "Generation", y = "Score",
             title = "Aneuploidy", col = "Sim label") +
        cinsim_theme()

      heterogeneity <- karyoSim$pop_measures_gw %>%
        filter(parameter == "heterogeneity") %>%
        ggplot(aes(x = g, group = sim_label, col = sim_label)) +
        geom_ribbon(aes(ymin = mean_value - sd_value,
                        ymax = mean_value + sd_value,
                        fill = sim_label),
                    alpha = 0.25, col = NA, show.legend = FALSE) +
        geom_line(aes(y = mean_value), size = 1) +
        labs(x = "Generation", y = "Score",
             title = "Heterogeneity", col = "Sim label") +
        cinsim_theme()

      deviation <- karyoSim$pop_measures_gw %>%
        filter(parameter == "deviation") %>%
        ggplot(aes(x = g, group = sim_label, col = sim_label)) +
        geom_ribbon(aes(ymin = mean_value - sd_value,
                        ymax = mean_value + sd_value,
                        fill = sim_label),
                    alpha = 0.25, col = NA, show.legend = FALSE) +
        geom_line(aes(y = mean_value), size = 1) +
        labs(x = "Generation", y = "Score",
             title = "Deviation", col = "Sim label") +
        cinsim_theme()

    }

  } else {

    # plot data per single replicate
    cell_count <- karyoSim$gen_measures %>%
      filter(parameter == "true_cell_count") %>%
      ggplot(aes(x = g, group = sim_label, col = sim_label)) +
      geom_line(aes(y = value), size = 1) +
      labs(x = "Generation", y = "Population size",
           title = "Growth", col = "Sim label") +
      scale_y_log10(breaks = scales::trans_breaks("log10", function(x) {10^x}),
                    labels = scales::trans_format("log10", scales::math_format(10^.x))) +
      cinsim_theme()

    fraction_surviving <- karyoSim$gen_measures %>%
      filter(parameter == "fraction_surviving") %>%
      ggplot(aes(x = g, group = sim_label, col = sim_label)) +
      geom_line(aes(y = value), size = 1) +
      labs(x = "Generation", y = "Fraction",
           title = "Survival", col = "Sim label") +
      cinsim_theme()

    mean_num_chr <- karyoSim$gen_measures %>%
      filter(parameter == "mean_num_chr") %>%
      ggplot(aes(x = g, group = sim_label, col = sim_label)) +
      geom_line(aes(y = value), size = 1) +
      labs(x = "Generation", y = "Number",
           title = "Chromosome count", col = "Sim label") +
      cinsim_theme()

    # plot data per chromosome
    if(per_chromosome) {

      aneuploidy <- karyoSim$pop_measures_chrom %>%
        filter(parameter == "aneuploidy") %>%
        ggplot(aes(x = g, group = chromosome, col = chromosome)) +
        geom_line(aes(y = value), size = 1) +
        facet_wrap(~sim_label) +
        scale_color_viridis_d() +
        scale_fill_viridis_d() +
        labs(x = "Generation", y = "Score",
             title = "Aneuploidy", col = "Chromosome") +
        cinsim_theme()

      heterogeneity <- karyoSim$pop_measures_chrom %>%
        filter(parameter == "heterogeneity") %>%
        ggplot(aes(x = g, group = chromosome, col = chromosome)) +
        geom_line(aes(y = value), size = 1) +
        facet_wrap(~sim_label) +
        scale_color_viridis_d() +
        scale_fill_viridis_d() +
        labs(x = "Generation", y = "Score",
             title = "Heterogeneity", col = "Chromosome") +
        cinsim_theme()

      deviation <- karyoSim$pop_measures_chrom %>%
        filter(parameter == "deviation") %>%
        ggplot(aes(x = g, group = chromosome, col = chromosome)) +
        geom_line(aes(y = value), size = 1) +
        facet_wrap(~sim_label) +
        scale_color_viridis_d() +
        scale_fill_viridis_d() +
        labs(x = "Generation", y = "Score",
             title = "Deviation", col = "Chromosome") +
        cinsim_theme()

      # plot data per genome-wide

    } else {

      aneuploidy <- karyoSim$pop_measures_gw %>%
        filter(parameter == "aneuploidy") %>%
        ggplot(aes(x = g, group = sim_label, col = sim_label)) +
        geom_line(aes(y = value), size = 1) +
        labs(x = "Generation", y = "Score",
             title = "Aneuploidy", col = "Sim label") +
        cinsim_theme()

      heterogeneity <- karyoSim$pop_measures_gw %>%
        filter(parameter == "heterogeneity") %>%
        ggplot(aes(x = g, group = sim_label, col = sim_label)) +
        geom_line(aes(y = value), size = 1) +
        labs(x = "Generation", y = "Score",
             title = "Heterogeneity", col = "Sim label") +
        cinsim_theme()

      deviation <- karyoSim$pop_measures_gw %>%
        filter(parameter == "deviation") %>%
        ggplot(aes(x = g, group = sim_label, col = sim_label)) +
        geom_line(aes(y = value), size = 1) +
        labs(x = "Generation", y = "Score",
             title = "Deviation", col = "Sim label") +
        cinsim_theme()

    }

  }

  # make plot list
  plot_list <- list(cell_count = cell_count,
                    fraction_surviving = fraction_surviving,
                    mean_num_chr = mean_num_chr,
                    aneuploidy = aneuploidy,
                    heterogeneity = heterogeneity,
                    deviation = deviation)

  if(plot) {
    p <- patchwork::wrap_plots(plot_list, nrow = 3, ncol = 2, byrow = FALSE) *
      theme(legend.position = "none")
  } else {
    return(plot_list)
  }

}
