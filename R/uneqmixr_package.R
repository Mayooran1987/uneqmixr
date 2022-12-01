#' Probability estimations and graphical displays in modelling the quantity of material sampled in the risk assessment.
#'
#'
#' @description This package aims to develop for getting probability estimations and graphical displays in the study associated with modelling the quantity of material sampled in the risk assessment.
#'
#' @details
#'
#'This package aims to develop probability estimations and graphical displays in modelling the quantity of material sampled in the risk assessment.
#'This study mainly focuses on the risk assessment when aggregating unequal incremental samples in the production process.
#'It mainly focuses on the risk assessment based on compound Poisson mixture distributions to model in five different scenarios.
#'
#'The following scenarios are considered for this study:
#'\enumerate{
#' \item Scenario 1— lots with homogeneous contamination;
#' \item Scenario 2— lots with heterogeneous, high-level contamination;
#' \item Scenario 3— lots with heterogeneous, low-level contamination;
#' \item Scenario 4 — lots with homogeneous contamination and concentration levels fluctuating from sub lots; and
#' \item Scenario 5—lots with heterogeneous contamination and concentration levels fluctuating from sub-lots.
#' }
#'
#' This package allows practitioners to get probability estimations and graphical displays based on this study. Also, this package can be used to validate the derived results in this study by simulation.
#'
#' @name uneqmixr
#'
#' @useDynLib uneqmixr, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#'
#' @importFrom stats rlnorm
NULL
