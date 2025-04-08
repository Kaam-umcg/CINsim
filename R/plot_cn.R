#' Plot CN frequencies over time
#'
#' This functions will plot the frequencie of copy numbers over time.
#'
#' @param karyoSim A karyoSim or karyoSimParallel object.
#' @param final_g A logical whether to plot the CN frequencies at the final generation.
#' @return A plot of copy number frequencies.
#' @author Bjorn Bakker
#' @import tidyverse
#' @export

plot_cn <- function(karyoSim, final_g = TRUE) {

  # check user input
  if(!any(class(karyoSim) %in% c("karyoSim", "karyoSimParallel", "karyoSim_simple", "karyoSimParallel_simple"))) {
    stop("An object of class karyoSim or karyoSimParallel is required.")
  }

  if(any(class(karyoSim) %in% c("karyoSim_simple", "karyoSimParallel_simple"))) {
    plot_cn_simple(karyoSim)
  } else {

    # compile tidysim if necessary
    if(class(karyoSim) %in% c("karyoSimParallel")) {
      karyoSim <- compile_karyoSimParallel(karyoSim, for_summary = FALSE)
      summarize <- TRUE
    } else {
      summarize <- FALSE
    }

    # get population measures
    pop_measures <- karyoSim$pop_measures

    if(final_g) {

      if(summarize) {

        # get cn frequencies at final generation
        cn <- pop_measures %>%
          group_by(sim_label) %>%
          filter(g == max(g)) %>%
          ungroup() %>%
          unnest(cn_freq)

        # summarize the data to get mean frequencies and sds
        cn <- cn %>%
          mutate(copy = factor(copy)) %>%
          group_by(chromosome, copy) %>%
          summarize(mean_freq = mean(freq, na.rm = TRUE),
                    sd_freq = sd(freq, na.rm = TRUE)) %>%
          ungroup() %>%
          replace_na(list(mean_freq = 0,
                          sd_freq = 0)) %>%
          group_by(chromosome) %>%
          mutate(mean_freq = mean_freq/sum(mean_freq),
                 cum_mean_freq = cumsum(mean_freq),
                 copy = factor(copy, levels = rev(sort(unique(cn$copy))))) %>%
          ungroup()

        # plot bar graphs with sds
        p <- cn %>%
          ggplot(aes(x = chromosome)) +
          geom_bar(aes(y = mean_freq, group = copy, fill = copy), stat = "identity") +
          geom_errorbar(aes(ymin = cum_mean_freq,
                            ymax = cum_mean_freq + sd_freq,
                            col = copy), width = 0.5) +
          scale_y_continuous(breaks = seq(0, 1, 0.2)) +
          scale_fill_manual(values = copy_num_cols) +
          scale_colour_manual(values = copy_num_cols, guide = FALSE) +
          coord_cartesian(ylim = c(0, 1)) +
          labs(x = "Chromosome", y = "Frequency", fill = "Copy") +
          cinsim_theme()

      } else {

        # get measures
        cn <- pop_measures %>%
          filter(g == max(g)) %>%
          ungroup() %>%
          unnest(cn_freq)

        # plot bar graphs wihout
        p <- cn %>%
          mutate(copy = factor(copy, levels = rev(sort(unique(cn$copy))))) %>%
          ggplot(aes(x = chromosome, y = freq)) +
          geom_bar(aes(fill = copy), stat = "identity") +
          scale_y_continuous(breaks = seq(0, 1, 0.2)) +
          scale_fill_manual(values = copy_num_cols) +
          coord_cartesian(ylim = c(0, 1)) +
          labs(x = "Chromosome", y = "Frequency", fill = "Copy") +
          cinsim_theme()

      }

      return(p)

    } else {

      # get cn frequencies at final generation
      cn <- pop_measures %>%
        unnest(cn_freq)

      if(summarize) {

        # summarize the data to get mean frequencies and sds
        cn <- cn %>%
          mutate(copy = factor(copy)) %>%
          group_by(g, chromosome, copy) %>%
          summarize(mean_freq = mean(freq, na.rm = TRUE),
                    sd_freq = sd(freq, na.rm = TRUE)) %>%
          ungroup() %>%
          replace_na(list(mean_freq = 0,
                          sd_freq = 0)) %>%
          group_by(g, chromosome) %>%
          mutate(mean_freq = mean_freq/sum(mean_freq),
                 cum_mean_freq = cumsum(mean_freq),
                 copy = factor(copy, levels = rev(sort(unique(cn$copy))))) %>%
          ungroup()

        # plot bar graphs with sds
        p <- cn %>%
          ggplot(aes(x = g)) +
          geom_bar(aes(y = mean_freq, group = copy, fill = copy),
                   stat = "identity", width = 0.75) +
          geom_errorbar(aes(ymin = cum_mean_freq,
                            ymax = cum_mean_freq + sd_freq,
                            col = copy), width = 0.25) +
          facet_wrap(~chromosome) +
          scale_y_continuous(breaks = seq(0, 1, 0.2)) +
          scale_fill_manual(values = copy_num_cols) +
          scale_colour_manual(values = copy_num_cols, guide = FALSE) +
          coord_cartesian(ylim = c(0, 1)) +
          labs(x = "Generation", y = "Frequency", fill = "Copy") +
          cinsim_theme()

      } else {

        # plot bar graphs wihout
        p <- cn %>%
          mutate(copy = factor(copy, levels = rev(sort(unique(cn$copy))))) %>%
          ggplot(aes(x = g, y = freq)) +
          geom_bar(aes(fill = copy), stat = "identity", width = 0.75) +
          facet_wrap(~chromosome) +
          scale_y_continuous(breaks = seq(0, 1, 0.2)) +
          scale_fill_manual(values = copy_num_cols) +
          coord_cartesian(ylim = c(0, 1)) +
          labs(x = "Generation", y = "Frequency", fill = "Copy") +
          cinsim_theme()

      }

      return(p)

    }

  }

}

