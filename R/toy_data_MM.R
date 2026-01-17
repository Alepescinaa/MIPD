#' Toy dataset of panel data generated under the Markovian assumption, recording disease status
#' for each individual at each observed visit.
#'
#' @format
#' A data frame with 13,528 rows and 10 columns:
#' \describe{
#'   \item{patient_id}{Unique identifier for each patient.}
#'   \item{cov1, cov2, cov3}{One continuous and two binary covariates.}
#'   \item{dead}{Binary indicator for whether the patient is deceased.}
#'   \item{death_time}{Time of death if occurred or censoring time otherwise.}
#'   \item{onset}{Binary indicator for disease onset.}
#'   \item{onset_age}{Age at disease onset or death_time otherwise.}
#'   \item{age}{Patient's age at the specific visit.}
#'   \item{visits}{Indicator of the current visit.}
#' }
#' @source Synthetic data.
"toy_data_MM"
