#' Prepare dataset to perform multiple imputation of the disease history
#'
#' @param data A `data.table` or `data.frame` that must contain the following columns:
#' - `patient_id`: Unique identifier for each patient (numeric).
#' - `dead`: Binary indicator (0/1) for whether the patient is dead.
#' - `death_time`: Time of death if it has occured or censoring time otherwise (numeric).
#' - `onset`: Binary indicator (0/1) for disease onset.
#' - `onset_age`: Age at disease onset if it has occured or death_time otherwise (numeric).
#' - `age`: Patient's current age at that specific visit (numeric).
#' - `visits`: Indicator of the current visit (numeric).
#' The `data.frame` could contain extra columns with covariates values.
#'
#' @returns A `list` containing two `data.frame`s:
#'
#' - The first `data.frame` retains all columns from the input data while adding `last_bfo`, which records
#'   the patient's age at their last visit before diagnosis (`onset == 1`). It includes only the first row
#'   per patient.
#' - The second `data.frame` is an ordered copy of the input data.
#'
#' Both `data.frame`s have updated, sequentially ordered `patient_id` values.
#'
#' @export
#'
#' @examples
#' data(toy_data_MM)
#' prepare_dataset(toy_data_MM)

#' @importFrom magrittr %>%

prepare_dataset <- function(data){

  n_pats <- length( unique(data$patient_id))
  scheme_visits <- data
  scheme_data <- data
  original_index <- unique(scheme_data$patient_id)
  scheme_data$last_bfo <- rep(NA, nrow(scheme_data))
  for (i in original_index){
    index <- which(scheme_data$onset[scheme_data$patient_id==i]==1)[1]
    scheme_data$last_bfo[scheme_data$patient_id==i] <- scheme_data$age[scheme_data$patient_id==i][index-1]
  }
  for (i in original_index){
    if(any(scheme_data$onset[scheme_data$patient_id==i]==1))
      scheme_data$onset[scheme_data$patient_id==i] <- 1
  }

  scheme_data <- scheme_data[order(scheme_data$patient_id),]

  row_id <- scheme_visits %>%
    dplyr::group_by(patient_id) %>%
    dplyr::summarise(nrows = dplyr::n())

  scheme_visits$patient_id <- rep(1:n_pats, times=as.numeric(row_id$nrows))
  scheme_data$patient_id <- rep(1:n_pats, times=as.numeric(row_id$nrows))

  scheme_data <- scheme_data %>%
    dplyr::group_by(patient_id) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()

  return(list(scheme_data=scheme_data, scheme_visits=scheme_visits ))

}
