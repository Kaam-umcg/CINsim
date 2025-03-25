#' Gather population measures
#'
#' This function will collect scores for aneuploidy and heterogeneity, and determine the
#' fraction of cells deviating from the modal copy number.
#'
#' @param karyotypes A matrix with cells in rows and chromosomes in colums. Rows and colums must be labelled appropriately.
#' @param euploid_ref The euploid copy number (user for aneuploidy score calculations - 2 by default).
#' @param g The generation identifier to add to the data frame.
#' @return A data frame containing all relevant population measures of the simulation.
#' @export
#' @author Bjorn Bakker

population_measures <- function(karyotypes = NULL, euploid_ref = 2, g = NULL) {

  # check user input
  if(is.null(karyotypes) | !is.matrix(karyotypes)) {
    stop("Karyotypes in matrix format should be provided")
  } else if(is.null(g)) {
    g <- "unspecified"
  }

  # gather measures and compile into tidy data frame,
  # and factor chromosome
  karyotype_measures <- tibble(chromosome = colnames(karyotypes),
                               aneuploidy = qAneuploidy(karyotypes, euploid_ref),
                               heterogeneity = qHeterogeneity(karyotypes),
                               deviation = dev_from_mode(karyotypes)) %>%
    mutate(chromosome = factor(chromosome, levels = colnames(karyotypes))) %>%
    nest(measures = c(chromosome, aneuploidy, heterogeneity, deviation))

  # transform karyotypes into data tibble
  karyotypes_tibble <- karyotypes %>%
    as_tibble() %>%
    mutate(cell_id = rownames(karyotypes))

  # calculate chromosomal copy number frequencies
  cn_freq <- karyotypes_tibble %>%
    gather(colnames(karyotypes), key = "chromosome", value = "copy") %>%
    mutate(chromosome = factor(chromosome, levels = colnames(karyotypes))) %>%
    group_by(chromosome, copy) %>%
    dplyr::count() %>%
    ungroup() %>% group_by(chromosome) %>%
    mutate(freq = n/sum(n)) %>%
    ungroup() %>%
    nest(cn_freq = c(n, freq))

  # determine clonal frequencies
  clon_freq <- karyotypes_tibble %>%
    group_by(cell_id) %>%
    dplyr::count() %>%
    ungroup() %>%
    mutate(freq = n/sum(n)) %>%
    nest(clone_freq = c(n, freq))

  # combine all data into a single nested data frame
  pop_measures <- bind_cols(karyotype_measures, cn_freq, clon_freq) %>%
    mutate(g = g)

  # return karyotype measures
  return(pop_measures)

}


