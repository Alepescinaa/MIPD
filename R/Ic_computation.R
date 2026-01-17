#' Compute Confidence Intervals using Rubin's Rule
#'
#' This function calculates confidence intervals for model parameters using Rubin's rule
#' for multiple imputations to take into account within and between variability.
#'
#' @param avg_parameters A `matrix` containing the average parameter estimates across imputations.
#' @param all_fits A `list` of fitted models from multiple imputed datasets.
#'
#' @returns A `list` with the same number of elements as the transitions of the multi-state model,
#'          each containing a matrix with the lower and upper 95% confidence interval bounds
#'          for the parameters associated with a specific transition.
#' @export
#'
#' @examples
#' IC <- Ic_computation(avg_parameters_MM, all_fits_MM)
#' print(IC)
#'
Ic_computation <- function(avg_parameters, all_fits){
   m <- length(all_fits)

  # Let's compute parameters covariance with Rubin rule
  param_matrix <- list()
  for (i in 1:m) {
    param_matrix[[i]] <- sapply(all_fits[[i]], stats::coefficients)
  }

  # within variance (mean of covariance matrix of each imputation)
  U_list <- list()
  for (i in 1:m) {
    U_list[[i]] <- lapply(all_fits[[i]], function(fit) stats::vcov(fit))
  }

  U_bar <- list()
  for (i in 1:3) {
    U_bar[[i]] <- matrix(0, nrow = nrow(U_list[[1]][[i]]), ncol = ncol(U_list[[1]][[i]]))

    for (j in 1:m) {
      U_bar[[i]] <- U_bar[[i]] + U_list[[j]][[i]]
    }
    U_bar[[i]] <- U_bar[[i]] / m
  }

  # between variance (variance of parameters estimate)
  param_matrix_mod1 <- list()
  param_matrix_mod2 <- list()
  param_matrix_mod3 <- list()
  for (i in 1:m) {
    param_matrix_mod1[[i]] <- param_matrix[[i]][, 1]
    param_matrix_mod2[[i]] <- param_matrix[[i]][, 2]
    param_matrix_mod3[[i]] <- param_matrix[[i]][, 3]
  }

  num_params <- ncol(avg_parameters)

  B1 <- matrix(0, nrow = num_params, ncol = num_params)
  for (j in 1:m) {
    B1 <- B1 + (param_matrix_mod1[[j]] - avg_parameters[1, ]) %*% t(param_matrix_mod1[[j]] - avg_parameters[1, ])
  }
  B1 <- B1 / (m - 1)

  B2 <- matrix(0, nrow = num_params, ncol = num_params)
  for (j in 1:m) {
    B2 <- B2 + (param_matrix_mod2[[j]] - avg_parameters[2, ]) %*% t(param_matrix_mod2[[j]] - avg_parameters[2, ])
  }
  B2 <- B2 / (m - 1)

  B3 <- matrix(0, nrow = num_params, ncol = num_params)
  for (j in 1:m) {
    B3 <- B3 + (param_matrix_mod3[[j]] - avg_parameters[3, ]) %*% t(param_matrix_mod3[[j]] - avg_parameters[3, ])
  }
  B3 <- B3 / (m - 1)

  # Rubin's rule
  T1 <- U_bar[[1]] + (1 + 1/m) * B1
  T2 <- U_bar[[2]] + (1 + 1/m) * B2
  T3 <- U_bar[[3]] + (1 + 1/m) * B3

  # confidence interval computation
  compute_CI <- function(R,theta_bar){
    SE <- sqrt(diag(R))
    alpha <- 0.05
    z_value <- stats::qnorm(1 - alpha / 2)
    lower_bound <- theta_bar - z_value * SE
    upper_bound <- theta_bar + z_value * SE

    conf_int <- cbind(lower_bound, upper_bound)
    colnames(conf_int) <- c("L95%", "U95%")
    return(conf_int)
}

  CI_mod1 <- compute_CI(T1, avg_parameters[1, ])
  CI_mod2 <- compute_CI(T2, avg_parameters[2, ])
  CI_mod3 <- compute_CI(T3, avg_parameters[3, ])

  return(list(CI_mod1, CI_mod2, CI_mod3))
}


