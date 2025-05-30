% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_coefficients.R, R/score2prob.R
\name{get_coefficients}
\alias{get_coefficients}
\title{Get coefficients to determine karyotype survival}
\usage{
get_coefficients(
  selection_metric = NULL,
  euploid_copy = 2,
  euploid_survival = 0.9,
  max_survival = 1
)

get_coefficients(
  selection_metric = NULL,
  euploid_copy = 2,
  euploid_survival = 0.9,
  max_survival = 1
)
}
\arguments{
\item{selection_metric}{The copy number frequency matrix.}

\item{euploid_copy}{The number of copies corresponding to the euploid state.}

\item{euploid_survival}{The probability of survival for a euploid karyotype.}

\item{max_survival}{The probability of survival for the most optimal karyotype.}

\item{selection_mode}{A selection method, either "cn_based" or "rel_copy".}
}
\value{
A named vector with slope, intercept and constant coefficients.

A named vector with slope and intercept coefficients.
}
\description{
This function will generate linear fit coefficients (slope, intercept, and constant) to convert
a karyotype fitness score into a survival probability

This function will generate linear fit coefficients (slope, intercept) to convert
a karyotype fitness score into a survival probability (cn_based only).
}
\author{
Bjorn Bakker
}
