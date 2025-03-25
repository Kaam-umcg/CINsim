#' Plot clonal frequencies (tidy)
#'
#' This function will allow you to plot the frequencies of the clones over time.
#'
#' @param karyoSim A karyoSim or karyoSimParallel object that contains a summary matrix for clonality frequency per generations.
#' @param line_graph A logical whether to plot frequencies as a line graph or a bar graph.
#' @return A table or plot of the clonal frequencies.
#' @import reshape2
#' @import tidyverse
#' @import ggplot2
#' @import tibble
#' @import tidyr
#' @import dplyr
#' @import RColorBrewer
#' @author Bjorn Bakker
#' @export

plot_clonality <- function(karyoSim, line_graph = TRUE) {

  # check user input
  if(class(karyoSim) != "karyoSim" & class(karyoSim) != "karyoSimParallel") {
    stop("An object of class karyoSim or karyoSimParallel is required")
  }

  if(class(karyoSim) == "karyoSim") {

    # get clonal frequencies
    clon_freq <- karyoSim$pop_measures %>%
      unnest(clone_freq) %>%
      mutate(cell_id = gsub(x = cell_id, pattern = "cell_", replacement = "") %>%
               as.numeric() %>%
               factor())

    if(line_graph) {

      # get top 1 clone at final
      top_clone <- clon_freq %>%
        filter(g == max(g),
               freq == max(freq)) %>%
        .$cell_id

      # plot line graph per clone
      p <- clon_freq %>%
        ggplot(aes(x = g, y = freq)) +
        geom_line(aes(group = cell_id, col = cell_id), linetype = 2, size = 0.5) +
        geom_line(data = clon_freq %>%
                    filter(cell_id == top_clone),
                  aes(group = cell_id, col = cell_id), size = 1) +
        ggrepel::geom_label_repel(data = clon_freq %>%
                                    filter(cell_id == top_clone,
                                           g == max(g)),
                                  aes(label = cell_id)) +
        scale_color_viridis_d() +
        labs(x = "Generation", y = "Frequency",
             col = "Clone", title = "Clonal abundance") +
        cinsim_theme()
      if(clon_freq$cell_id %>% unique() %>% length() > 20) {
        p <- p + theme(legend.position = "none")
      }

    } else {

      # plot bar graph
      p <- clon_freq %>%
        ggplot(aes(x = g, y = freq)) +
        geom_bar(aes(fill = cell_id), col = "grey", stat = "identity") +
        scale_y_continuous(breaks = seq(0, 1, 0.2)) +
        scale_fill_viridis_d() +
        labs(x = "Generation", y = "Frequency",
             fill = "Clone", title = "Clonal abundance") +
        cinsim_theme()
      if(clon_freq$cell_id %>% unique() %>% length() > 20) {
        p <- p + theme(legend.position = "none")
      }

    }

    return(p)

  } else if(class(karyoSim) == "karyoSimParallel") {

    # get clonal frequencies
    clon_freq <- karyoSim %>%
      purrr::map_df("pop_measures") %>%
      unnest(clon_freq) %>%
      mutate(cell_id = gsub(x = cell_id, pattern = "cell_", replacement = "") %>%
               as.numeric())

    if(line_graph) {

      # get cummulative frequencies add rank and factor accordingly
      top_clones <- clon_freq %>%
        filter(g == max(g)) %>%
        group_by(cell_id) %>%
        summarize(median_freq = median(freq, na.rm = TRUE)) %>%
        mutate(rank = rank(-median_freq, ties.method = "min")) %>%
        arrange(rank) %>%
        top_n(n = -5, wt = rank) %>%
        mutate(cell_id = factor(cell_id, levels = .$cell_id))

      # summarize data
      clon_freq <- clon_freq %>%
        group_by(g, cell_id) %>%
        summarize(median_freq = median(freq, na.rm = TRUE),
                  min_freq = min(freq, na.rm = TRUE),
                  max_freq = max(freq, na.rm = TRUE)) %>%
        ungroup()

      top_clones_stat <- clon_freq %>%
        filter(cell_id %in% top_clones$cell_id) %>%
        mutate(cell_id = factor(cell_id, levels = levels(top_clones$cell_id)))


      # plot line graph per clone
      p <- clon_freq %>%
        filter(!cell_id %in% top_clones$cell_id) %>%
        ggplot(aes(x = g, y = median_freq, group = cell_id)) +
        #geom_linerange(aes(ymin = min_freq, ymax = max_freq), col = "grey") +
        geom_line(col = "grey", linetype = 2, size = 0.5, alpha = 0.5) +
        #geom_linerange(data = top_clones_stat,
        #               aes(ymin = min_freq, ymax = max_freq, col = cell_id)) +
        geom_line(data = top_clones_stat, aes(col = cell_id), size = 1) +
        ggrepel::geom_label_repel(data = top_clones_stat %>%
                                    filter(g == max(g)),
                                  aes(label = cell_id)) +
        scale_color_viridis_d(direction = -1) +
        labs(x = "Generation", y = "Median frequency",
             col = "Top 5 clones", title = "Clonal abundance") +
        cinsim_theme()

    } else {

      # summarize stats
      clon_freq_summary <- clon_freq %>%
        group_by(g, cell_id) %>%
        summarize(mean_freq = mean(freq, na.rm = TRUE),
                  sd_freq = sd(freq, na.rm = TRUE)) %>%
        ungroup() %>%
        replace_na(list(mean_freq = 0,
                        sd_freq = 0)) %>%
        group_by(g) %>%
        mutate(mean_freq = mean_freq/sum(mean_freq),
               cum_mean_freq = cumsum(mean_freq),
               cell_id = factor(cell_id, levels = rev(sort(unique(clon_freq$cell_id))))) %>%
        ungroup()

      # plot bar graph
      p <- clon_freq_summary %>%
        ggplot(aes(x = g, y = mean_freq)) +
        geom_bar(aes(fill = cell_id), col = "grey", stat = "identity") +
        geom_errorbar(aes(ymin = cum_mean_freq,
                          ymax = cum_mean_freq + sd_freq,
                          col = cell_id), width = 0.25) +
        scale_y_continuous(breaks = seq(0, 1, 0.2)) +
        scale_fill_viridis_d() +
        scale_colour_viridis_d(guide = FALSE) +
        labs(x = "Generation", y = "Frequency",
             fill = "Clone", title = "Clonal abundance") +
        cinsim_theme()

      if(clon_freq_summary$cell_id %>% unique() %>% length() > 20) {
        p <- p + theme(legend.position = "none")
      }

    }

    return(p)

  }

}
