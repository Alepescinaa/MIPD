fit_model_cox <- function (data_imputed,cov_vector)
{

  # Define transition matrix and the number of transitions for multi-state model
  tmat <- mstate::transMat(x = list(c(2, 3), c(3), c()), names = c("Disease-free", "Disease", "Death"))
  n_trans <- max(tmat, na.rm = TRUE)  # Get the total number of transitions

  covariate_names <- cov_vector
  # Prepare the data for the multi-state model using msprep function
  data_long <- mstate::msprep(data = data_imputed, trans = tmat,
                              time = c(NA, "onset_age", "death_time"),
                              status = c(NA, "onset", "dead"),
                              keep = c("age",covariate_names),
                              id = "patient_id")

  # Adjust the start times for transitions
  data_long$Tstart[data_long$trans < 3] <- data_long$Tstart[data_long$trans < 3] + data_long$age[data_long$trans < 3]
  data_long$time <- data_long$Tstop - data_long$Tstart

  # Expand covariates for the model
  data_long <- mstate::expand.covs(data_long, covariate_names)

  # Define covariates for Cox Proportional Hazards model
  expanded_covariates <- setdiff(names(data_long), c("patient_id", "from", "to", "trans", "Tstart", "Tstop", "time", "status", "age", covariate_names))

  # Create formula for Cox Proportional Hazards model
  formula_str <- paste("survival::Surv(Tstart, Tstop ,status) ~",
                       paste(expanded_covariates, collapse = " + "),
                       "+ strata(trans)")
  model_formula <- stats::as.formula(formula_str)

  # Fit Cox Proportional Hazards model using Breslow method
  model_cox <- survival::coxph(model_formula, data = data_long, method = "breslow")



  # Identify coefficient names and their indices
  coef_names <- names(model_cox$coefficients)


  # Identify coefficients for the different transitions
  to_onset_indices <- grep("\\.1$", coef_names)
  to_death_indices <- grep("\\.2$", coef_names)
  dis_to_death_indices <- grep("\\.3$", coef_names)

  # Extract corresponding coefficients for each transition
  onsetpar <- model_cox$coefficients[to_onset_indices]
  deathpar <- model_cox$coefficients[to_death_indices]
  deathpar_dis <- model_cox$coefficients[dis_to_death_indices]



  return(list(onsetpar, deathpar, deathpar_dis))

}
