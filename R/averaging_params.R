#' Averaging parameter estimates across all fitted models for each transition
#'
#' This function averages the parameter estimates obtained from multiple fitted models
#' across each transition. It returns a `data.frame` where each row corresponds to
#' a transition and each column represents a parameter, containing the mean value
#' of each parameter across all provided models.
#'
#' @param all_fits A `list` of previously fitted models, where each element in the list
#' is itself a list containing the model coefficients for each transition.
#'
#' @returns A `data.frame` with rows representing transitions and columns representing
#' parameters. The values in the `data.frame` are the averaged parameter estimates
#' across different models.
#'
#' @export
#'
#' @examples
#' avg_params <- averaging_params(all_fits_MM)
#' print(avg_params)

averaging_params <- function(all_fits) {

  # Extract parameter names from the first fitted model
  param_names <- names(stats::coef(all_fits[[1]][[1]]))
  n_params <- length(param_names)

  # Initialize matrix for averaged parameters
  averaged_params <- matrix(0, nrow = length(all_fits[[1]]), ncol = n_params)
  colnames(averaged_params) <- param_names

  # Loop over transitions and parameters
  for (i in seq_along(all_fits[[1]])) {
    for (p in seq_along(param_names)) {
      averaged_params[i, p] <- mean(sapply(all_fits, function(fit) stats::coef(fit[[i]])[p]))
    }
  }

  return(averaged_params)
}
