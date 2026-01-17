library(purrr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

size <- 5000

ci_root <- file.path(
  "~/OneDrive - Politecnico di Milano/ANNO 5/PAPER/MIPD/CI",
  paste0("high_prev_", size)
)


scheme_ids <- c(1,2,3,4,5,6)
term_rate <- "high"
models <- c ("cox", "flexsurv_panel", "flexsurv_real", "msm", "msm_age", "nhm", "mipd")
params <- c ("cov")

out_root <- file.path(
  "~/OneDrive - Politecnico di Milano/ANNO 5/PAPER/MIPD/CI",
  paste0("high_prev_", size)
)


all_ci <- map_dfr(scheme_ids, function(sid) {
  base <- file.path(ci_root, sprintf("cov_scheme_%d", sid), term_rate)
  subdirs <- list.dirs(base, full.names = TRUE, recursive = FALSE)

  map_dfr(subdirs, function(sd) {
    fp <- file.path(sd, "ic_summary.rds")
    if (!file.exists(fp)) return(tibble())

    readRDS(fp) %>%
      mutate(
        scheme_id = sid
      ) %>%

      mutate(
        est_mean = round(est_mean,3),
        # rel_bias
        coverage = scales::percent(coverage_rate, accuracy = 1),
        coverage_rate = NULL,
        type_I   = scales::percent(type_I_rate,   accuracy = 1),
        type_I_rate = NULL,
        type_II  = scales::percent(type_II_rate,  accuracy = 1),
        type_II_rate = NULL
      ) %>%
      filter(model %in% models, param %in% params)
  })
})


table <- function(pm_data){
  table <- split(pm_data, pm_data$trans_idx)
  table <- lapply(table, function(x)
    x %>% dplyr::select(scheme_id, param, model, est_mean, coverage_rate, type_I_rate, type_II_rate)) #add rel_bias

  table <- lapply(table, function(x)
    x %>% mutate_at(3:5, round, 3) %>% mutate_at(6, round, 2)  %>%
      mutate_all(as.character) %>%
      filter(parameter%in%c("beta_cov1","beta_cov2","beta_cov3")))

  table <- do.call("rbind",lapply(1:3,function(y){
    res <- split(table[[y]],table[[y]]$parameter)
    tab <- do.call("rbind",lapply(1:3, function(x){
      res[[x]] %<>% dplyr::select(-parameter)
      res2 <- rbind(rep(paste0("**beta_",x,"**"),5),res[[x]] %>% slice(c(5,2,1,4,3)))
      return(res2)
    }))
    tab <- cbind(tab,trans=rep(y,nrow(tab)))
  }

  ))

}
