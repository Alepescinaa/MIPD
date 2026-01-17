#' @noRd
.pool_rubin <- function(Q, U_list) {
  # Q: m x p matrix of estimates
  # U_list: list length m of p x p covariance matrices
  m <- nrow(Q)
  Qbar <- colMeans(Q)

  Ubar <- Reduce(`+`, U_list) / m
  B <- stats::cov(Q)
  Tcov <- Ubar + (1 + 1/m) * B       # Rubin total variance

  # Barnard-Rubin df per-parameter
  diagU <- diag(Ubar)
  diagB <- diag(B)
  denom <- (1 + 1/m) * diagB
  df <- rep(Inf, length(Qbar))
  ok <- denom > 0
  df[ok] <- (m - 1) * (1 + diagU[ok] / denom[ok])^2

  list(Qbar = Qbar, Ubar = Ubar, B = B, Tcov = Tcov, df = df, m = m)
}

.strip_flexsurvreg <- function(obj, keep = c(
  "call", "ncovs", "ncoveffs", "basepars", "covpars", "concat.formula",
  "cov", "coefficients", "npars", "fixedpars", "optpars"
)) {
  for (nm in setdiff(names(obj), keep)) {
    obj[[nm]] <- NULL
  }
  obj
}

.extract_metadata <- function(fit) {
  list(
    distribution     = distribution,
    clock_assumption = clock_assumption,
    cov_vector       = cov_vector,
    custom_formula   = custom_formula
  )
}

pool_rubin_one_model <- function(fits, cl = 0.95) {

  fits <- Filter(function(f) inherits(f, "flexsurvreg"), fits)
  if (length(fits) < 2)
    stop("Need at least 2 successful fits to pool.")

  par_names <- names(stats::coef(fits[[1]]))
  if (any(vapply(
    fits,
    function(f) !identical(names(stats::coef(f)), par_names),
    logical(1)
  ))) {
    stop("Parameter names/order differ across imputations.")
  }

  Q <- do.call(rbind, lapply(fits, stats::coef))
  U_list <- lapply(fits, stats::vcov)

  rub <- .pool_rubin(Q, U_list)

  pooled <- .strip_flexsurvreg(fits[[1]])
  pooled$coefficients <- rub$Qbar
  pooled$cov          <- rub$Tcov

  attr(pooled, "rubin") <- rub
  class(pooled) <- "flexsurvreg_pooled"


  pooled
}

pool_rubin_all_transitions <- function(
    all_fits, cl = 0.95,
    distribution, clock_assumption, cov_vector, custom_formula,
    loss_plots = NULL,
    logs_cols  = NULL
) {


  metadata <- list(
    distribution     = distribution,
    clock_assumption = clock_assumption,
    cov_vector       = cov_vector,
    custom_formula   = custom_formula
  )

  # optional extras
  if (!is.null(loss_plots)) metadata$loss_plots <- loss_plots
  if (!is.null(logs_cols))  metadata$logs_cols  <- logs_cols

  ok <- which(vapply(all_fits, Negate(is.null), logical(1)))
  if (!length(ok)) stop("No successful fits found.")
  template <- all_fits[[ok[1]]]

  # --- single-transition fit
  if (inherits(template, "flexsurvreg")) {
    pooled <- pool_rubin_one_model(all_fits, cl = cl)
    attr(pooled, "metadata") <- metadata
    return(pooled)
  }

  # --- multi-state: list of flexsurvreg
  if (is.list(template) && all(vapply(template, inherits, logical(1), "flexsurvreg"))) {

    K <- length(template)
    pooled_list <- vector("list", K)

    for (k in seq_len(K)) {
      kth_fits <- lapply(all_fits, function(obj) if (is.list(obj)) obj[[k]] else NULL)
      pooled_list[[k]] <- pool_rubin_one_model(kth_fits, cl = cl)
    }

    names(pooled_list) <- names(template)
    class(pooled_list) <- "flexsurvreg_pooled_multistate"

    # attach ALL metadata (including loss_plots if present)
    attr(pooled_list, "metadata") <- metadata
    return(pooled_list)
  }

  stop("Unsupported fit structure.")
}
