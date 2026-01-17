#' @noRd
.get_rubin <- function(x) {
  rub <- attr(x, "rubin", exact = TRUE)
  if (is.null(rub)) stop("Missing Rubin diagnostics in attr(x, 'rubin').")
  rub
}


#------------------------------------------------------------
# coef / vcov
#------------------------------------------------------------
#' @noRd
#' @exportS3Method stats::coef flexsurvreg_pooled
coef.flexsurvreg_pooled <- function(object, ...) {
  .get_rubin(object)$Qbar
}
#' @noRd
#' @exportS3Method stats::vcov flexsurvreg_pooled
vcov.flexsurvreg_pooled <- function(object, ...) {
  .get_rubin(object)$Tcov
}

#------------------------------------------------------------
# confint with Rubin df (t-interval per-parameter)
#------------------------------------------------------------

#' @noRd
#' @exportS3Method stats::confint flexsurvreg_pooled
confint.flexsurvreg_pooled <- function(object, parm = NULL, level = 0.95, ...) {

  rub <- .get_rubin(object)

  est <- rub$Qbar
  V   <- rub$Tcov
  se  <- sqrt(diag(V))
  df  <- rub$df
  names(df) <- names(est)

  if (!is.null(parm)) {
    est <- est[parm]
    se  <- se[parm]
    df  <- df[parm]
  }

  alpha <- 1 - level
  crit <- vapply(df, function(d) {
    if (is.finite(d)) stats::qt(1 - alpha/2, df = d)
    else stats::qnorm(1 - alpha/2)
  }, numeric(1))

  out <- cbind(
    lower = est - crit * se,
    upper = est + crit * se
  )

  attr(out, "level") <- level
  out
}

#------------------------------------------------------------
# print
#------------------------------------------------------------

#' @noRd
#' @exportS3Method base::print flexsurvreg_pooled
print.flexsurvreg_pooled <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  rub <- .get_rubin(x)
  est <- rub$Qbar
  se  <- sqrt(diag(rub$Tcov))
  df  <- rub$df
  level <- 0.95

  ci <- confint.flexsurvreg_pooled(x, level = level)

  tab <- data.frame(
    est = est,
    se  = se,
    df  = df,
    LCL = ci[,1],
    UCL = ci[,2],
    row.names = names(est)
  )

  cat("Pooled flexsurvreg (Rubin's rules)\n")
  if (!is.null(x$dist)) cat("Distribution:", x$dist, "\n")
  if (!is.null(rub$m))   cat("m imputations:", rub$m, "\n\n")

  print(round(tab, digits = digits))
  invisible(x)
}

#------------------------------------------------------------
# summary flexsurvreg_pooled:
# Rubin table (t-based CI)
#------------------------------------------------------------
#' @noRd
#' @exportS3Method stats::summary flexsurvreg_pooled
summary.flexsurvreg_pooled <- function(object, level = 0.95, ...) {

  rub <- .get_rubin(object)

  est <- rub$Qbar
  se  <- sqrt(diag(rub$Tcov))
  df  <- rub$df

  alpha <- 1 - level
  crit <- vapply(df, function(d) {
    if (is.finite(d)) stats::qt(1 - alpha/2, df = d)
    else stats::qnorm(1 - alpha/2)
  }, numeric(1))

  tab <- data.frame(
    est = est,
    se  = se,
    df  = df,
    LCL = est - crit * se,
    UCL = est + crit * se,
    row.names = names(est)
  )

  out <- list(
    call  = object$call,
    table = tab,
    rubin = rub,
    level = level
  )

  class(out) <- "summary.flexsurvreg_pooled"
  out
}

#------------------------------------------------------------
# summary flexsurvreg_pooled_multistate:
#------------------------------------------------------------
#' @noRd
#' @exportS3Method stats::summary flexsurvreg_pooled_multistate
summary.flexsurvreg_pooled_multistate <- function(
    object,
    digits = max(3L, getOption("digits") - 3L),
    plots_loss = FALSE,
    plots = plots_loss,
    ...
) {
  metadata <- attr(object, "metadata", exact = TRUE)

  m <- NULL
  if (length(object) > 0) {
    rub <- attr(object[[1]], "rubin", exact = TRUE)
    if (!is.null(rub) && !is.null(rub$m)) m <- rub$m
  }

  cat("\nPooled multi-state model (flexsurv)\n")
  cat(strrep("-", nchar("Pooled multi-state model (flexsurv)")), "\n", sep = "")
  cat("Each transition is fitted as a `flexsurvreg_pooled` model; parameters are pooled using Rubin's rules.\n")
  if (!is.null(m)) cat(sprintf("Imputations (m): %d\n", m))
  cat(sprintf("Transitions: %d\n", length(object)))

  # ---- metadata
  if (!is.null(metadata) && length(metadata) > 0) {
    cat("\nMetadata\n")
    cat(strrep("-", nchar("Metadata")), "\n", sep = "")

    # print everything except loss_plots
    printed_any <- FALSE
    for (nm in names(metadata)) {
      if (nm == "loss_plots") next
      printed_any <- TRUE
      cat(sprintf("  %-20s: %s\n", nm, paste(deparse(metadata[[nm]]), collapse = " ")))
    }
    if (!printed_any) cat("  (no scalar metadata fields; only `loss_plots` stored)\n")
  }

  # ---- optional plots
  show_plots <- isTRUE(plots_loss) || isTRUE(plots)
  if (show_plots && !is.null(metadata) && !is.null(metadata$loss_plots)) {
    cat("\nLoss diagnostics plots\n")
    cat(strrep("-", nchar("Loss diagnostics plots")), "\n", sep = "")

    lp <- metadata$loss_plots
    if (!is.null(lp$raw$total)) print(lp$raw$total)
    if (!is.null(lp$raw$aligned_plus_fixmatch)) print(lp$raw$aligned_plus_fixmatch)
    if (!is.null(lp$weighted$total)) print(lp$weighted$total)
    if (!is.null(lp$weighted$aligned_plus_fixmatch)) print(lp$weighted$aligned_plus_fixmatch)
  }

  cat("\n")
  invisible(object)
}
