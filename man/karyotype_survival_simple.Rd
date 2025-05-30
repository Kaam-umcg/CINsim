% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Cinsim_simple.R
\name{karyotype_survival_simple}
\alias{karyotype_survival_simple}
\title{Determine cell survival (simplified)}
\usage{
karyotype_survival_simple(
  karyotypes,
  copy_num_boundaries = c(1, 8),
  selection_mode = NULL,
  selection_metric = NULL,
  a = NULL,
  b = NULL,
  c = NULL,
  qMod = NULL
)
}
\arguments{
\item{karyotypes}{A matrix with cells in rows and chromosomes in colums. Rows and colums must be labelled appropriately.}

\item{copy_num_boundaries}{A vector of length 2 with the minimum and maximum viable copy number state for any chromosome set.}

\item{selection_mode}{One of three possible modes, "cn_based", "rel_copy", "davoli", for copy number-based, relative copy numbed-based or Davoli-based selection.}

\item{selection_metric}{A fitness metric for the cn_based, rel_copy or davoli selection measure.}

\item{a}{Value of exponent for karyotype selection.}

\item{b}{Value of constant for karyotype selection.}

\item{c}{Value of second constant for karyotype selection.}

\item{qMod}{A modifier for the survival probability.}
}
\value{
A vector of logicals whether the cells have survived or not.
}
\description{
This function will determine the survival probability for a karyotype, depending on specific parameters.
}
\author{
Bjorn Bakker
}
