#' Multiple imputation for interval-censored illness-death model and fitting Parametric Illness-Death Model over imputed datasets
#'
#' Performs multiple imputation on a dataset containing interval-censored data
#' about the disease progression of individuals, and fits a parametric illness-death model
#' assuming exact observation times over each dataset with imputed disease status and disease onset.
#'
#' @param data A `data.table` or `data.frame` containing the following columns:
#' - `patient_id`: Unique identifier for each patient (numeric).
#' - `dead`: Binary indicator (0/1) of whether the patient is dead.
#' - `death_time`: Time of death if it has occurred, or censoring time otherwise (numeric).
#' - `onset`: Binary indicator (0/1) of disease onset.
#' - `onset_age`: Age at disease onset if it has occurred, or death time otherwise (numeric).
#' - `age`: Patient's current age at the specific visit (numeric).
#' - `visits`: Indicator of the current visit (numeric).
#' Additional columns in the `data.frame` will be treated as covariates.
#'
#' @param cov_vector vector of covariates names.
#' @param m A numeric value specifying the number of imputations to perform.
#' We recommend choosing a value between 20 and 50 to ensure consistency.
#' Keep in mind that the computational time increases with the number of imputations (`m`),
#' so for large datasets, it may be more efficient to choose a smaller value like `m = 30`.
#'
#' @param clock_assumption A character string specifying the assumed time scale:
#'   - `"forward"` (Markov assumption): Only the time at transition is considered.
#'   - `"mix"` (Semi-Markov assumption): Accounts for the time spent in the disease-free state.
#'
#' @param distribution A character string specifying the parametric form of the baseline hazard.
#'   It must be one of the distributions available in `flexsurv::flexsurvreg`.
#'
#' @param max_iter maximum number of iterations.
#' @param eps tolerance for convergence.
#' @returns A `list` containing the average values of model parameter estimates in the first element
#'  and 95% confidence intervals computed using Rubin's rule saved by transition in the second element.
#'
#' @export
#'
#' @examples
#' # Example usage:
#' fit <- run_mipd(
#'toy_data_MM,
#'m = 5,
#'cov_vector = c("cov1", "cov2", "cov3"),
#'clock_assumption = "forward",
#'distribution = "gompertz",
#'inner_cores = 1,
#'max_iter = 20,
#'eps = 1e-3
#')
#' )
#' print(fit)

run_mipd <- function(data, m, clock_assumption, cov_vector, distribution,  inner_cores, max_iter, eps, boot = FALSE){

  # ---- Only change from your original: mclapply wrapper (Windows-safe) ----

  .mc_apply <- function(X, FUN, ...) {
    if (.Platform$OS.type == "windows") {
      # mclapply uses forking and is not available on Windows; fall back
      lapply(X, FUN, ...)
    } else {
      parallel::mclapply(X, FUN, mc.cores = inner_cores, ...)
    }
  }
  # -------------------------------------------------------------------------

  check_dataset_structure <- function(data) {

    if (!inherits(data, c("data.table", "data.frame"))) {
      stop("The input data must be a data.table or data.frame.")
    }

    required_columns <- c("patient_id", "dead", "death_time", "onset", "onset_age", "age", "visits",cov_vector)

    missing_columns <- setdiff(required_columns, colnames(data))
    if (length(missing_columns) > 0) {
      stop("The following required columns are missing: ", paste(missing_columns, collapse = ", "))
    }

    if (!is.numeric(data$patient_id)) {
      stop("The 'patient_id' column must be numeric.")
    }
    if (!all(data$dead %in% c(0, 1))) {
      stop("The 'dead' column must be binary.")
    }
    if (!is.numeric(data$death_time)) {
      stop("The 'death_time' column must be numeric.")
    }
    if (!all(data$onset %in% c(0, 1))) {
      stop("The 'onset' column must be binary.")
    }
    if (!is.numeric(data$onset_age)) {
      stop("The 'onset_age' column must be numeric.")
    }
    if (!is.numeric(data$age)) {
      stop("The 'age' column must be numeric.")
    }
    if (!is.numeric(data$visits)) {
      stop("The 'visits' column must be numeric.")
    }
    na_columns <- names(data)[sapply(data, anyNA)]
    if (length(na_columns) > 0) {
      stop("The following columns contain missing values (NA): ", paste(na_columns, collapse = ", "))
    }

    return(TRUE)
  }

  check_clock_assumption <- function(clock_assumption) {

    if (!is.character(clock_assumption) || length(clock_assumption) != 1) {
      stop("The 'clock_assumption' parameter must be a single character string.")
    }

    if (!clock_assumption %in% c("forward", "mix")) {
      stop("The 'clock_assumption' parameter must be either 'forward' or 'mix'.")
    }

    return(TRUE)
  }

  check_distribution <- function(distribution) {

    if (!is.character(distribution) || length(distribution) != 1) {
      stop("The 'distribution' parameter must be a single character string.")
    }

    supported_distributions <- c("genf", "genf.orig", "gengamma", "gengamma.orig",
                                 "exp", "weibull", "weibullph", "lnorm", "gamma",
                                 "gompertz", "llogis", "exponential", "lognormal")

    if (!(distribution %in% supported_distributions)) {
      stop(paste("The 'distribution' parameter must be one of the following:",
                 paste(supported_distributions, collapse = ", ")))
    }

    return(TRUE)
  }

  check_dataset_structure(data)
  check_clock_assumption(clock_assumption)
  check_distribution(distribution)

  # Step 1: Create base visit and per-patient data
  n_pats <- length(unique(data$patient_id))
  scheme_visits <- data %>% dplyr::arrange(patient_id, visits)  # keep order stable
  scheme_data   <- scheme_visits
  original_index <- unique(scheme_data$patient_id)

  # ensure numeric vector
  scheme_data$last_bfo <- rep(NA_real_, nrow(scheme_data))

  for (i in original_index) {
    onset_idx <- which(scheme_data$onset[scheme_data$patient_id == i] == 1)[1]
    if (!is.na(onset_idx) && onset_idx > 1) {
      scheme_data$last_bfo[scheme_data$patient_id == i] <-
        scheme_data$age[scheme_data$patient_id == i][onset_idx - 1]
    }
  }

  # Mark patient as having onset = 1 if any visit has onset
  for (i in original_index) {
    if (any(scheme_data$onset[scheme_data$patient_id == i] == 1)) {
      scheme_data$onset[scheme_data$patient_id == i] <- 1
    }
  }

  scheme_data <- scheme_data[order(scheme_data$patient_id), ]

  # keep row counts aligned with patient ordering
  row_id <- scheme_visits %>%
    dplyr::count(patient_id, name = "nrows") %>%
    dplyr::arrange(patient_id)

  scheme_visits$patient_id <- rep(seq_len(n_pats), times = as.numeric(row_id$nrows))
  scheme_data$patient_id   <- rep(seq_len(n_pats), times = as.numeric(row_id$nrows))

  scheme_data <- scheme_data %>%
    dplyr::group_by(patient_id) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::select(-visits)

  # Step 2: Define last_bfo correctly
  last_visits <- scheme_visits %>%
    dplyr::group_by(patient_id) %>%
    dplyr::slice_max(visits, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(patient_id, visits, last_visit_age = age)

  scheme_data <- merge(scheme_data, last_visits, by = "patient_id", all.x = TRUE, sort = FALSE)
  scheme_data$last_bfo <- ifelse(scheme_data$onset == 1, scheme_data$last_bfo, scheme_data$last_visit_age)
  scheme_data <- scheme_data %>% dplyr::select(-last_visit_age)

  # Step 3: Define covariates and prepare output
  scheme_data <- as.data.frame(scheme_data)
  covariate_names <- cov_vector


  temp <- scheme_data

  # Define transition matrix and the number of transitions for multi-state model
  tmat <- mstate::transMat(x = list(c(2, 3), c(3), c()), names = c("Disease-free", "Disease", "Death"))
  n_trans <- max(tmat, na.rm = TRUE)  # Get the total number of transitions


  # Prepare the data for the multi-state model using msprep function
  data_long <- mstate::msprep(data = temp, trans = tmat,
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


  # Compute baseline hazards for each transition
  all_haz <- survival::basehaz(model_cox, center = FALSE)

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

  # Extract hazard values for each transition
  onsethaz <- all_haz[all_haz$strata == "trans=1", 1:2]
  deathhaz <- all_haz[all_haz$strata == "trans=2", 1:2]
  deathhaz_dis <- all_haz[all_haz$strata == "trans=3", 1:2]


  #Compute hazard
  compute_hazard <- function(df) {
    df <- rbind(data.frame(time = 0, hazard = 0), df)
    df$h <- c(diff(df$hazard) / diff(df$time), NA)
    return(df)
  }

  # Apply to all three hazard tables
  onsethaz <- compute_hazard(onsethaz)
  deathhaz <- compute_hazard(deathhaz)
  deathhaz_dis <- compute_hazard(deathhaz_dis)

  # Merge onset and death hazards
  combhaz <- merge(onsethaz, deathhaz, by = "time", all = TRUE, suffixes = c(".12", ".13"))

  # Ensure first values are not NA
  combhaz$hazard.12[1] <- ifelse(is.na(combhaz$hazard.12[1]), 0, combhaz$hazard.12[1])
  combhaz$hazard.13[1] <- ifelse(is.na(combhaz$hazard.13[1]), 0, combhaz$hazard.13[1])

  # Precompute time differences
  dt <- c(NA, diff(combhaz$time))

  # Fill NAs in hazard rates and recalculate cumulative hazards
  na_12 <- is.na(combhaz$h.12)
  na_13 <- is.na(combhaz$h.13)

  for (i in which(na_12)) {
    combhaz$h.12[i] <- combhaz$h.12[i - 1]
    combhaz$hazard.12[i] <- combhaz$hazard.12[i - 1] + combhaz$h.12[i - 1] * dt[i]
  }

  for (i in which(na_13)) {
    combhaz$h.13[i] <- combhaz$h.13[i - 1]
    combhaz$hazard.13[i] <- combhaz$hazard.13[i - 1] + combhaz$h.13[i - 1] * dt[i]
  }

  # Matrices to store imputed data
  disease_status <- matrix(rep(scheme_data$onset, m), nrow = nrow(scheme_data), ncol = m)
  disease_age <- matrix(0, nrow = nrow(scheme_data), ncol = m)
  set.seed(2)
  c1 <- matrix(stats::runif(nrow(scheme_data) * m), nrow = nrow(scheme_data), ncol = m)
  c2 <- matrix(stats::runif(nrow(scheme_data) * m), nrow = nrow(scheme_data), ncol = m)

  # Covariate matrix
  covar <- as.matrix(scheme_data[, covariate_names])

  # Linear predictor calculations
  lp.12 <- covar %*% as.vector(t(onsetpar))
  lp.13 <- covar %*% as.vector(t(deathpar))
  lp.23 <- covar %*% as.vector(t(deathpar_dis))

    criteria <- Inf
    old <- list(onsetpar, deathpar, deathpar_dis)
    k <- 0
    while (criteria > eps && k < max_iter) {


      n_patients <- nrow(scheme_data)

      # Windows‑compatible parallel backend
      future::plan(future::multisession, workers = inner_cores)
      results <- future.apply::future_lapply(
        X = 1:n_patients,
        FUN = function(x) {
          MIPD:::process_patient(
            x,
            combhaz$time,
            combhaz$hazard.12,
            combhaz$hazard.13,
            combhaz$h.12,
            lp.12,
            lp.13,
            c1,
            c2,
            scheme_visits[scheme_visits$patient_id == x, ],
            m,
            mean(scheme_data$onset)
          )
        },future.seed=T
      )

      disease_age <- do.call(rbind, lapply(results, `[[`, "age"))
      disease_status <- do.call(rbind, lapply(results, `[[`, "status"))


      result_iteration <- parallel::mclapply(1:m, function(j) {
        temp <- scheme_data
        disease_generated <- which(disease_status[, j] == 1)
        temp$onset[disease_generated] <- 1
        temp$onset_age[disease_generated] <- disease_age[disease_generated, j]
        fit_model_cox(temp, cov_vector)
      })


      average_coefficients <- function(lst) {
        n <- 3
        lapply(seq_len(n), function(i) {
          components <- lapply(lst, `[[`, i)
          Reduce("+", components) / length(components)
        })
      }

      avg_coefficients <- average_coefficients(result_iteration)

      criteria1 <- norm(as.matrix((avg_coefficients[[1]] - old[[1]]) / old[[1]]))
      criteria2 <- norm(as.matrix((avg_coefficients[[2]] - old[[2]]) / old[[2]]))
      criteria3 <- norm(as.matrix((avg_coefficients[[2]] - old[[2]]) / old[[2]]))
      criteria <- criteria1+criteria2+criteria3


      lp.12 <- covar %*% avg_coefficients[[1]]
      lp.13 <- covar %*% avg_coefficients[[2]]
      lp.23 <- covar %*% avg_coefficients[[3]]

      print(avg_coefficients)
      old <- avg_coefficients

      k <- k+1
      message(sprintf(
        "Iteration %d — Total criteria: %.5f (c1: %.5f, c2: %.5f,  c3: %.5f)",
        k, criteria, criteria1, criteria2, criteria3
      ))
    }


    all_fits <-  parallel::mclapply(1:m, function(j) {
      temp <- scheme_data
      disease_generated <- which(disease_status[, j] == 1)
      temp$onset[disease_generated] <- 1
      temp$onset_age[disease_generated] <- disease_age[disease_generated, j]
      fit_model(temp, cov_vector, clock_assumption, distribution)
    })



    # Compute averaged parameters
    #averaged_params <- averaging_params(all_fits)

    # Bootstrap

    # bias_hat = NULL
    # B <- 100
    # m_boot <- 5
    #
    # if(boot == TRUE){
    #   boots <- bias_bootstrap(
    #     scheme_data = scheme_data,
    #     disease_status = disease_status,
    #     disease_age = disease_age,
    #     cov_vector = cov_vector,
    #     clock_assumption = clock_assumption,
    #     distribution = distribution,
    #     averaging_params = averaging_params,
    #     B = B, m_boot = m_boot, seed = 42
    #   )
    #
    #   boot_mean <- Reduce("+", boots) / length(boots)
    #   bias_hat <- boot_mean - averaged_params
    # }


    #-------------------------------------
    # Pooling with Rubin's rules
    #-------------------------------------

    pooled_fit <- pool_rubin_all_transitions(all_fits, cl = 0.95, distribution, clock_assumption, cov_vector, custom_formula = NULL)

    return(pooled_fit)
}
