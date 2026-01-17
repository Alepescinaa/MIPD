#' Fits a parametric illness-death model with known disease progression
#'
#' This function fits a parametric illness-death model leveraging the `mstate` and `flexsurv` packages,
#' assuming either a Markov or Semi-Markov structure for disease progression.
#'
#' @param data A `data.frame` containing imputed disease onset times, disease status, and covariates for each individual.
#' @param cov_vector vector with covariates names contained in data.
#' @param clock_assumption A character string specifying the assumed time scale:
#'   - `"forward"` (Markov assumption): Only the time at transition is considered.
#'   - `"mix"` (Semi-Markov assumption): Accounts for the time spent in the disease-free state.
#' @param distribution A character string specifying the parametric form of baseline hazards.
#'   Must be one of the distributions available in `flexsurv::flexsurvreg`, e.g., `"weibull"`, `"exponential"`, `"gompertz"`.
#'
#' @returns A `list` of fitted `flexsurvreg` models, one per transition in the illness-death model.
#'
#' @export
#'
#' @examples
#' # Example usage:
#' data <- imputed_dataset_MM # Example dataset with required variables
#' fit <- fit_model(data,c("cov1","cov2","cov3"), clock_assumption = "forward", distribution = "weibull")
#' summary(fit[[1]]) # Summary of first transition model

fit_model <- function(data, cov_vector, clock_assumption, distribution) {
  tmat <- mstate::transMat(x = list(c(2, 3), c(3), c()), names = c("Disease-free", "Disease", "Death"))


  covariate_names <- cov_vector


  data_long <- mstate::msprep(data = data, trans = tmat,
                              time = c(NA, "onset_age", "death_time"),
                              status = c(NA, "onset", "dead"),
                              keep = c("age", covariate_names),
                              id = "patient_id")

  data_long$Tstart[data_long$trans < 3] <- data_long$Tstart[data_long$trans < 3] + data_long$age[data_long$trans < 3]
  data_long$time <- data_long$Tstop - data_long$Tstart

  n_trans <- max(tmat, na.rm = TRUE)
  fits <- vector(mode = "list", length = n_trans)

  forward_formula <- stats::as.formula(paste("survival::Surv(Tstart, Tstop, status) ~", paste(covariate_names, collapse = " + ")))
  reset_formula <- stats::as.formula(paste("survival::Surv(time, status) ~", paste(covariate_names, collapse = " + ")))


  if (clock_assumption == "forward") {
    for (i in 1:3) {
      fits[[i]] <- flexsurv::flexsurvreg(forward_formula,
                                         data = subset(data_long, trans == i),
                                         dist = distribution)
    }
  } else if (clock_assumption == "mix") {
    for (i in 1:2) {
      fits[[i]] <- flexsurv::flexsurvreg(forward_formula,
                                         data = subset(data_long, trans == i),
                                         dist = distribution)
    }
    fits[[3]] <- flexsurv::flexsurvreg(reset_formula,
                                       data = subset(data_long, trans == 3),
                                       dist = distribution)
  }

  return(fits)
}
