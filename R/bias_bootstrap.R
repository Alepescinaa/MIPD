bias_bootstrap <- function(
    scheme_data, disease_status, disease_age,
    cov_vector, clock_assumption, distribution,
    averaging_params,
    B = 100, m_boot = 5L, seed = 1L
) {
  set.seed(seed)

  n_patients <- nrow(scheme_data)
  m <- ncol(disease_status)
  stopifnot(n_patients == nrow(disease_age), m == ncol(disease_age))

  use_j <- sort(sample.int(m, size = min(m_boot, m), replace = FALSE))

  boots <- vector("list", B)

  for (b in seq_len(B)) {

    samp_idx <- sample.int(n_patients, size = n_patients, replace = TRUE)

    fits_b <- lapply(use_j, function(j) {
      temp_b <- scheme_data[samp_idx, , drop = FALSE]

      temp_b$onset     <- disease_status[samp_idx, j]
      temp_b$onset_age <- disease_age[samp_idx, j]


      if ("patient_id" %in% names(temp_b)) {
        temp_b$patient_id <- seq_len(nrow(temp_b))
      }

      fit_model(temp_b, cov_vector, clock_assumption, distribution)
    })

    # Pool across the m_boot imputations for this bootstrap draw → 3 x p pooled estimate
    boots[[b]] <- averaging_params(fits_b)
  }

  boots
}
