#' Plot a copy number heatmap.
#'
#' This function will plot a copy number heatmap based on a karyotype matrix.
#'
#' @param karyoSim A karyoSim object (final output of CINsim).
#' @param subset_size The number of cells to subset and plot from the karyotype matrix.
#' @param clones A character vector of the clone identifiers to plot.
#' @param file A (optional) name for the pdf-file to which the plot will be printed.
#' @return A heatmap with copy numbers shown in colours.
#' @import tidyverse
#' @import ggplot2
#' @author Bjorn Bakker
#' @export

cnvHeatmap <- function(karyoSim = NULL, subset_size = 1000, clones = NULL, file = NULL) {
  # check user input
  if (!class(karyoSim) %in% c("karyoSim", "karyoSim_simple")) {
    message("An object of class karyoSim is required")
    stop
  }

  # subselect karyotype data frame and get chromosome names
  kldf_plot <- karyoSim$karyotypes
  chromosomes <- colnames(kldf_plot)

  # add cell and clone identifiers
  kldf_plot <- kldf_plot %>%
    as.data.frame() %>%
    tibble::rowid_to_column("cell_id") %>%
    mutate(clone_id = rownames(kldf_plot) %>%
      gsub(x = ., pattern = "cell_", replacement = "") %>%
      as.numeric())

  # subset the tibble for a particular clone if it is provided
  if (is.null(clones)) {
    clones <- unique(kldf_plot$clone_id)
    subsetted <- FALSE
  } else if (!(clones %in% unique(kldf_plot$clone_id))) {
    message("Clone ID not present in karyotype matrix - plotting all clones")
    clones <- unique(kldf_plot$clone_id)
    subsetted <- FALSE
  } else {
    subsetted <- TRUE
  }

  kldf_plot <- kldf_plot %>%
    filter(clone_id %in% clones)

  # randomly sample a subset when number of cells exceeds the subset size
  if (nrow(kldf_plot) > subset_size) {
    kldf_plot <- kldf_plot %>%
      filter(cell_id %in% sample(kldf_plot$cell_id, size = subset_size, replace = FALSE))
  }

  # reorder the karyotypes based on similarity
  ord <- kldf_plot %>%
    select(chromosomes) %>%
    dist(method = "euclidean") %>%
    hclust(method = "ward.D") %>%
    .$order

  # melt the cnv table and make a tibble; also factor the cell identifier based on similarity
  kldf_plot <- kldf_plot %>%
    gather(chromosomes, key = "chromosome", value = "copy") %>%
    mutate(
      cell_id = factor(cell_id, levels = kldf_plot$cell_id[ord]),
      chromosome = factor(chromosome, levels = chromosomes)
    ) %>%
    as_tibble()

  # set plot title
  plot_title <- paste0(
    "pMisseg: ", karyoSim$sim_info["pMisseg"], ", generations: ", karyoSim$sim_info[["g"]],
    ", number of cells: ", length(unique(kldf_plot$cell_id)), "/", nrow(karyoSim$karyotypes)
  )
  if (subsetted) {
    plot_title <- paste(plot_title, ",\n subsetted for clone(s) ", paste0(clones, sep = " "), sep = "")
  }

  # plot the cnv heatmap
  g <- kldf_plot %>%
    ggplot(aes(x = chromosome, y = cell_id)) +
    geom_tile(aes(fill = factor(copy))) +
    scale_fill_manual(values = copy_num_cols) +
    labs(x = "Chromosomes", y = "Cells", fill = "Copy", title = plot_title) +
    guides(fill = guide_legend(ncol = 2)) +
    cinsim_theme() +
    theme(
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank()
    )

  # save to file if requested, otherwise rerturn the plot
  if (is.null(file)) {
    return(g)
  } else {
    ggsave(plot = g, filename = paste0(file, ".pdf"), width = 10, height = 8)
  }
}
