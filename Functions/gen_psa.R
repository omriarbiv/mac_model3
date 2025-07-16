library(dampack)
library(tidyverse)
source("microsim.R")

gen_psa_params <- function(n_psa = 1000, seed = 1, pr = 0) {
  
  ### Baseline values
  # Obtaining baseline values
  baseline.vals <- read.csv("values/baseline.csv", 
                            row.names = "variable_name")
  params_base <- as.list(baseline.vals$value)
  params_name <- rownames(baseline.vals)
  names(params_base) <- params_name
  
  life.table <- read.csv("values/2022_lifetable.csv") |>
    select("Age", "qx") |>  # qx is the death rate
    rename(age = Age, base_rate = qx) |>
    # Changing age to be just the number
    mutate(age = as.numeric(str_split_i(age, " ", i = 1)))
   
  # # Retrieves the rate of death from the life table (bcc = base culture 
  # # conversion)
  # params_base$p_death <- life.table |> 
  #   # Adding additional rows so that we can align ages
  #   add_row(age = 111:120, base_rate = rep(1, 10)) |> 
  #   # Death rate for those who had culture conversion in first 3 months
  #   mutate(cc6m_rate = base_rate * params_base$hr_cc6m)
  
  psa_params <- params_base
  
  #############################################################################
  ### PSA parameters
  # number of iterations
  n_psa <- 1000 
  # prior (i.e., Jeffreys prior)
  pr <- 0 
  set.seed(2)
  
  # Utilities
  # From Shah et al 2021
  u1 <- dampack::beta_params(0.82, 0.17)
  psa_params$u_obv <- rbeta(n_psa, u1$alpha, u1$beta)
  psa_params$u_T <- rbeta(n_psa, u1$alpha, u1$beta)
  
  # From Ahmadian et al 2022
  psa_params$DisU_SE_Azi <- rbeta(n_psa, 66343, 2475535)
  # Assuming similar disutility as azithromycin
  psa_params$DisU_SE_Rif <- rbeta(n_psa, 66343, 2475535)
  
  # From Law et al 2014 and Brown et al 2001
  u2 <- dampack::beta_params(0.4, 0.29)
  psa_params$U_SE_Etb <- rbeta(n_psa, u2$alpha, u2$beta)
  
  # From Park et al 2021
  u3 <- dampack::beta_params(0.93, 0.1902)
  psa_params$u_cure <- rbeta(n_psa, u3$alpha, u3$beta)
  
  ## Probabilities
  # From Jeong et al 2015
  psa_params$p_trtd1 <- r_p(p_r(rbeta(n_psa, 7 + pr, 118 - 7 + pr)), 6)
  psa_params$p_trtd2 <- r_p(p_r(rbeta(n_psa, 3 + pr, 118 - 3 + pr)), 6)
  psa_params$p_vis <- r_p(p_r(rbeta(n_psa, 1 + pr, 118 - 1 + pr)), 24)
  
  # From Wallace et al 2014
  psa_params$p_cure_6mo <- r_p(p_r(rbeta(n_psa, 154 + pr, 180 - 154 + pr), 6))
  # From CONVERT trial
  psa_params$p_cure_late <- r_p(p_r(rbeta(n_psa, 65 + pr, 224 - 65 + pr), 6))
  
  # From Kwon et al 20202
  psa_params$p_cure_wo_emb_6mo <- r_p(p_r(rbeta(n_psa, 32 + pr, 47 - 32 + pr), 6))
  tmp1_rr <- r_p(32/47, 6) / r_p(154 / 180, 6)
  psa_params$p_cure_wo_emb_late <- r_p(p_r(rbeta(n_psa,  65 + pr, 
                                                 224 - 65 + pr), 6) * tmp1_rr)
  # Same as ethambutol, assumption
  psa_params$p_cure_wo_rif_6mo <- r_p(p_r(rbeta(n_psa, 32 + pr, 47 - 32 + pr), 6))
  psa_params$p_cure_wo_rif_late <- r_p(p_r(rbeta(n_psa, 65 + pr, 
                                                 224 - 65 + pr), 6) * tmp1_rr)
  
  # From Kadota et al 2016
  psa_params$p_cure_wo_azi_6mo <- r_p(p_r(rbeta(n_psa, 12 + pr, 33 - 12 + pr), 12))
  psa_params$p_cure_wo_azi_late <- r_p(p_r(rbeta(n_psa, 12 + pr, 33 - 12 + pr), 12))
  
  
  # Selecting names of probabilities/utilities except death
  # From Moon et al 2019
  psa_params$p_spont_cure <- r_p(p_r(rbeta(n_psa, 26 + pr, 157 - 26 + pr), 43.2))
  # From Hwang et al 2017
  psa_params$p_prog <- r_p(p_r(rbeta(n_psa, 22 + pr, 115 - 22 + pr), 31.2))
  # From Wallace et al 2014
  psa_params$p_recurrence <- r_p(p_r(rbeta(n_psa, 65 + pr, 402 - 65 + pr), 42))
  
  # Hazard ratio
  # Unpublished data
  hrv <- ((baseline.vals["hr_cc6m", "max_value"] - 
             baseline.vals["hr_cc6m", "min_value"]) / 3.92) ^ 2
  hrp <- lnorm_params(params_base$hr_cc6m, hrv)
  psa_params$hr_cc6m <- rlnorm(n_psa, meanlog = hrp$mu, sdlog = hrp$sigma)
  
  # unchanged variables, not in PSA
  psa_params$start_age <- rep(params_base$start_age, n_psa)
  psa_params$discount <- rep(params_base$discount, n_psa)
  psa_params$u_death <- rep(params_base$u_death, n_psa)
  psa_params$rif_discount <- NULL
  psa_params$emb_discount <- NULL
  
  # Initializing results variable
  psa_results <- data.frame(
    Run = NA,
    Strategies = NA,
    QALYs = NA
  )
  psa_names <- names(psa_params)
  psa_det_obs <- list()
  psa_det_treat <- list()
  
}
