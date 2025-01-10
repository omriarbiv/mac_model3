library(tidyverse)
# library(dampack)
# library(darthtools)

###########################################################
### Helper functions
r_p <- function(rate, time = 1) { return(1 - exp(-rate * time)) }  
p_r <- function(prob, time = 1) { return(-log(1 - prob) / time) } 
source("Functions/probs.R")
source("Functions/effs.R")

###########################################################
### Microsim function
microsim_model <- function(l_params = params_base, 
                           arm = "treat", 
                           n_p = 1000, fup = 40, 
                           cycle_length = 1 / 12) {
  
  # set.seed(given_seed)
  n_c <- fup / cycle_length     
  # All states
  v_s <- c("Obs",       # Observation
           "T",         # First line treatment
           "TMV",       # Modify treatment for EMB due to visual S/E 
           "TM1",       # Modify treatment for RMP due to S/E
           "TM2",       # Modify treatment for AZI due to S/E
           "TM1+2",     # Modify treatment for both
           "TMV+TM1",   # Modify treatment for EMB and RMP due to S/E
           "TMV+TM2",   # Modify treatment for EMB and AZI due to S/E
           "TMV+TM1+2", # Modify treatment for EMB, RMP, and AZI 
           "Cure",      # Cure
           "Death")
  
  # Counter variables
  ## prev_cc6m = prior culture conversion in 6 months
  ## prev_vis_se = prior visual side effects
  ## obsT counter = number of times gone from obs to treatment
  ## lit = length in treatment
  prev_cc6m <- prev_vis_se <- rep(FALSE, n_p)
  obsTcounter <- lit <- empty_vec <- rep(0, n_p)
  empty_tm <- matrix(data = 0, nrow = n_p, ncol = length(v_s),
                     dimnames = list("patient" = 1:n_p, "state" = v_s))
  state_matrix <- util_matrix <- matrix(
    data = NA, nrow = n_p, ncol = n_c + 1,
    dimnames = list("patient" = 1:n_p, "cycle" = 0:n_c))
  
  if(arm == "treat") {
    state_matrix[, 1] <- "T"
    util_matrix[, 1] <- l_params$u_T
  } else if (arm == "observe") {
    state_matrix[, 1] <- "Obs"
    util_matrix[, 1] <- l_params$u_obv
  } else {
    stop("Inadequate 'arm'. Either 'treat' or 'observe'")
  }
  
  for (cy in 1:n_c) {
    
    probs_mat <- probs(n_p = n_p, cycle = cy,
                      state_matrix, v_s, l_params, lit,
                      prev_cc6m, prev_vis_se, empty_tm,
                      empty_vec)
    
    state_matrix[, cy + 1] <- sapply(1:n_p, \(x) 
                                     sample(v_s, 
                                            size = 1, 
                                            prob = probs_mat[x, ]))
    # state_matrix[, cy + 1] <- samplev(probs_mat, 1) --> bug with the samplev function!
    
    tmp <- effs(state_matrix, cy, l_params,
                prev_cc6m, prev_vis_se, obsTcounter, lit,
                empty_vec)
    util_matrix[, cy + 1] <- tmp$util
    prev_cc6m <- tmp$prev_cc6m
    prev_vis_se <- tmp$prev_vis_se
    obsTcounter <- tmp$obsTcounter
    lit <- tmp$lit
    
  }
  
  ## Applying discounting
  discount_vec <- 1 / ((1 + (l_params$discount / 12)) ^ (0:n_c))
  util_total <- util_matrix %*% discount_vec
  
  return(
    list(
      "util_total" = util_total,
      "state_matrix" = state_matrix,
      "util_matrix" = util_matrix,
      "prev_cc6m" = prev_cc6m,
      "prev_vis_se" = prev_vis_se,
      "obsTcounter" = obsTcounter
    )
  )
}
