library(dampack)
library(tidyverse)
source("microsim.R")

gen_psa_params <- function(n_psa, seed = 1, pr = 0,
                           baseline.vals = baseline.vals, 
                           psa_params = params_base, 
                           life.table = life.table) {
  
  
  ### PSA parameters
  # For debugging
  # n_psa <- 10
  # psa_params <- params_base
  set.seed(seed)
  
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
  
  # From Park et al 2021, after meta-analysis
  u3 <- dampack::beta_params(0.904, 0.016)
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
  # From Miwa et al
  psa_params$p_cure_wo_rif_6mo <- r_p(p_r(rbeta(n_psa, 29 + pr, 40-29 + pr), 6))
  psa_params$p_cure_wo_rif_late <- r_p(p_r(rbeta(n_psa, 4 + pr, 
                                                 40 - 4 + pr), 6))
  
  # Continue on azitrho
  psa_params$p_cure_wo_azi_6mo <- psa_params$p_cure_6mo
  psa_params$p_cure_wo_azi_late <- psa_params$p_cure_late
  
  
  # Selecting names of probabilities/utilities except death
  # From Moon et al 2019
  psa_params$p_spont_cure <- r_p(p_r(
    rbeta(n_psa, 26 + pr, 157 - 26 + pr), 
    43.2))
  # From Kwon et al
  psa_params$p_prog <- r_p(p_r(
    rbeta(n_psa, 25 + pr, 228 - 25 + pr), 
    33.6))
  # From Wallace et al 2014
  psa_params$p_recurrence <- r_p(p_r(
    rbeta(n_psa, 65 + pr, 402 - 65 + pr), 42))
  
  # Hazard ratio
  hr_cc6m <- rlnorm(n_psa, 
                    meanlog = psa_params$hrm1, 
                    sdlog = sqrt(psa_params$hrv1))
  hr <- rlnorm(n_psa, 
               meanlog = psa_params$hrm2, 
               sdlog = sqrt(psa_params$hrv2))
  
  # psa_params$p_death <- replicate(n_psa, list(life.table))
  ntm_rate <- lapply(hr, \(x) x * life.table[,2])
  cc6m_rate <- lapply(1:n_psa, \(x) ntm_rate[[x]] * hr_cc6m[[x]])
  psa_params$p_death <- lapply(1:n_psa, \(x) 
                               cbind(life.table, 
                                     ntm_rate = ntm_rate[[x]], 
                                     cc6m_rate = cc6m_rate[[x]]))
  
  # unchanged variables, not in PSA
  psa_params$start_age <- rep(psa_params$start_age, n_psa)
  psa_params$discount <- rep(psa_params$discount, n_psa)
  psa_params$u_death <- rep(psa_params$u_death, n_psa)
  
  # New vars
  psa_params$DisU_T <- rep(0, n_psa)
  
  # From Hong et al 2014
  u5 <- dampack::beta_params((0.966-0.949), sqrt(0.015^2 + 0.017^2))
  psa_params$DisU_obsT <- rbeta(n_psa, u5$alpha, u5$beta)
  
  return(psa_params)
  
}
