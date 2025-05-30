% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/score2prob.R
\name{make_cinsim_coeffcients}
\alias{make_cinsim_coeffcients}
\title{Get list of coefficients to determine karyotype survival}
\usage{
make_cinsim_coeffcients(
  selection_metric = NULL,
  euploid_copy = 2,
  min_survival_euploid = 0.9,
  max_survival_euploid = 0.9,
  max_survival = 1,
  interval = 0.1,
  probability_types = c("pSurvival")
)
}
\arguments{
\item{selection_metric}{The copy number frequency matrix.}

\item{euploid_copy}{The number of copies corresponding to the euploid state.}

\item{min_survival_euploid}{The minimum probability of survival for a euploid karyotype.}

\item{max_survival_euploid}{The maximum probability of survival for a euploid karyotype.}

\item{max_survival}{The maximum probability of survival for the optimal karyotype.}

\item{interval}{The interval of the survival probabilities (range min to max).}

\item{probability_types}{The probabilities over which to iterate (either one of pDivision, pMisseg, pSurvival)}
}
\value{
A list of lists with with slope and intercept coefficients.
}
\description{
This function will generate linear fit coefficients (slope, intercept) to convert
a karyotype fitness score into a survival probability (cn_based only).
}
\author{
Bjorn Bakker
}
