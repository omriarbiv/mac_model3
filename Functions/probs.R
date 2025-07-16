probs <- function(n_p, cycle, state_matrix, v_s,
                  r, lit, prev_cc6m, prev_vis_se, tm,
                  p_death2) {
  
  ## r = l_params
  ## p_death2 = death_rate
  
  # Choose the right death rate and divide it into monthly at a
  # constant rate based on current age
  cur_age <- round(r$start_age) + floor(cycle / 12)
  
  # Getting monthly rate for given age
  dage <- r$p_death[r$p_death[, "age"] == cur_age] / 12
  # Changing to probability depending on if prior culture conversion in 6m
  p_death2[prev_cc6m == T] <- r_p(dage[4])
  p_death2[prev_cc6m == F] <- r_p(dage[3])
  
  cs <- state_matrix[, cycle]
  
  ## Death
  tm[cs == "Death", "Death"] <- 1
  
  ## Treatment
  ##############################################################################
  # First 6 months
  ##############################################################################
  ## On standard treatment
  ### Probability of rif intolerance
  tm[lit <= 6 & cs == "T", "TM1"] <- 
    (1 - p_death2[cs == "T" & lit <= 6]) * 
    r$p_trtd1 *
    (1 - r$p_trtd2) * (1 - r$p_vis)
  ### Probability of azi intolerance
  tm[lit <= 6 & cs == "T", "TM2"] <- 
    (1 - p_death2[cs == "T" & lit <= 6]) * 
    r$p_trtd2 *
    (1 - r$p_trtd1) * (1 - r$p_vis)
  ### Probability of emb intolerance
  tm[lit <= 6 & cs == "T", "TMV"] <- 
    (1 - p_death2[cs == "T" & lit <= 6]) * 
    r$p_vis *
    (1 - r$p_trtd1) * (1 - r$p_trtd2)
  ### Probability of death
  tm[lit <= 6 & cs == "T", "Death"] <- 
    p_death2[cs == "T" & lit <= 6] * 
    (1 - r$p_vis) *
    (1 - r$p_trtd1) * (1 - r$p_trtd2)
  ### Probability of remaining on treatment
  tm[lit <= 6 & cs == "T", "T"] <- 
    (1 - p_death2[cs == "T" & lit <= 6]) * 
    (1 - r$p_vis) *
    (1 - r$p_trtd1) * (1 - r$p_trtd2)
  
  ## On TM1
  tm[lit <= 6 & cs == "TM1", "TM1+2"] <- 
    r$p_trtd2 * 
    (1 - r$p_vis) * 
    (1 - p_death2[cs == "TM1" & lit <= 6])
  tm[lit <= 6 & cs == "TM1", "TMV+TM1"] <- 
    r$p_vis * 
    (1 - r$p_trtd2) * 
    (1 - p_death2[cs == "TM1" & lit <= 6])
  tm[lit <= 6 & cs == "TM1", "Death"] <- 
    p_death2[cs == "TM1" & lit <= 6] * 
    (1 - r$p_trtd2) * 
    (1 - r$p_vis)
  tm[lit <= 6 & cs == "TM1", "TM1"] <- 
    (1 - p_death2[cs == "TM1" & lit <= 6]) * 
    (1 - r$p_trtd2) * 
    (1 - r$p_vis)
  
  ## On TM2
  tm[lit <= 6 & cs == "TM2", "TM1+2"] <- 
    r$p_trtd1 * 
    (1 - r$p_vis) * 
    (1 - p_death2[cs == "TM2" & lit <= 6])
  tm[lit <= 6 & cs == "TM2", "TMV+TM2"] <- 
    r$p_vis * 
    (1 - r$p_trtd1) * 
    (1 - p_death2[cs == "TM2" & lit <= 6])
  tm[lit <= 6 & cs == "TM2", "Death"] <- 
    p_death2[cs == "TM2" & lit <= 6] * 
    (1 - r$p_trtd1) * 
    (1 - r$p_vis)
  tm[lit <= 6 & cs == "TM2", "TM2"] <- 
    (1 - p_death2[cs == "TM2" & lit <= 6]) * 
    (1 - r$p_trtd1) * 
    (1 - r$p_vis)
  
  ## On TMV
  tm[lit <= 6 & cs == "TMV", "TMV+TM1"] <- 
    r$p_trtd1 * 
    (1 - r$p_trtd2) * 
    (1 - p_death2[cs == "TMV" & lit <= 6])
  tm[lit <= 6 & cs == "TMV", "TMV+TM2"] <- 
    r$p_trtd2 * 
    (1 - r$p_trtd1) * 
    (1 - p_death2[cs == "TMV" & lit <= 6])
  tm[lit <= 6 & cs == "TMV", "Death"] <- 
    p_death2[cs == "TMV" & lit <= 6] * 
    (1 - r$p_trtd1) * 
    (1 - r$p_trtd2)
  tm[lit <= 6 & cs == "TMV", "TMV"] <- 
    (1 - p_death2[cs == "TMV" & lit <= 6]) * 
    (1 - r$p_trtd1) * 
    (1 - r$p_trtd2)
  
  ## On TM1+2
  tm[lit <= 6 & cs == "TM1+2", "TMV+TM1+2"] <- 
    r$p_vis * 
    (1 - p_death2[cs == "TM1+2" & lit <= 6])
  tm[lit <= 6 & cs == "TM1+2", "Death"] <- 
    p_death2[cs == "TM1+2" & lit <= 6] * 
    (1 - r$p_vis)
  tm[lit <= 6 & cs == "TM1+2", "TM1+2"] <- 
    (1 - p_death2[cs == "TM1+2" & lit <= 6]) * 
    (1 - r$p_vis)
  
  ## On TMV+TM1
  tm[lit <= 6 & cs == "TMV+TM1", "TMV+TM1+2"] <- 
    r$p_trtd2 * 
    (1 - p_death2[cs == "TMV+TM1" & lit <= 6])
  tm[lit <= 6 & cs == "TMV+TM1", "Death"] <- 
    p_death2[cs == "TMV+TM1" & lit <= 6] * 
    (1 - r$p_trtd2)
  tm[lit <= 6 & cs == "TMV+TM1", "TMV+TM1"] <- 
    (1 - p_death2[cs == "TMV+TM1" & lit <= 6]) * 
    (1 - r$p_trtd2)
  
  ## On TMV+TM2
  tm[lit <= 6 & cs == "TMV+TM2", "TMV+TM1+2"] <- 
    r$p_trtd1 * 
    (1 - p_death2[cs == "TMV+TM2" & lit <= 6])
  tm[lit <= 6 & cs == "TMV+TM2", "Death"] <- 
    p_death2[cs == "TMV+TM2" & lit <= 6] * 
    (1 - r$p_trtd1)
  tm[lit <= 6 & cs == "TMV+TM2", "TMV+TM2"] <-  
    (1 - p_death2[cs == "TMV+TM2" & lit <= 6]) * 
    (1 - r$p_trtd1) 
  
  ## On TMV+TM1+2
  ### Probability of death
  tm[lit <= 6 & cs == "TMV+TM1+2", "Death"] <- 
    p_death2[cs == "TMV+TM1+2" & lit <= 6]
  ### Probability of remaining on treatment
  tm[lit <= 6 & cs == "TMV+TM1+2", "TMV+TM1+2"] <- 
    (1 - p_death2[cs == "TMV+TM1+2" & lit <= 6])
  
  # Month 7-12
  ## On standard treatment
  ### Probability of emb intolerance
  tm[lit > 6 & lit < 13 & cs == "T", "TMV"] <- 
    r$p_vis * 
    (1 - p_death2[cs == "T" & lit > 6 & lit < 13])
  ### Probability of death
  tm[lit > 6 & lit < 13 & cs == "T", "Death"] <- 
    p_death2[cs == "T" & lit > 6 & lit < 13] * 
    (1 - r$p_vis)
  ### Probability of remaining on treatment
  tm[lit > 6 & lit < 13 & cs == "T", "T"] <- 
    (1 - p_death2[cs == "T" & lit > 6 & lit < 13]) * 
    (1 - r$p_vis)
  
  ## TM1
  ### Probability of emb intolerance
  tm[lit > 6 & lit < 13 & cs == "TM1", "TMV+TM1"] <- 
    r$p_vis * 
    (1 - p_death2[cs == "TM1" & lit > 6 & lit < 13])
  ### Probability of death
  tm[lit > 6 & lit < 13 & cs == "TM1", "Death"] <- 
    p_death2[cs == "TM1" & lit > 6 & lit < 13] * 
    (1 - r$p_vis)
  ### Probability of remaining on treatment
  tm[lit > 6 & lit < 13 & cs == "TM1", "TM1"] <-  
    (1 - p_death2[cs == "TM1" & lit > 6 & lit < 13]) * 
    (1 - r$p_vis)
  
  ## TM2
  ### Probability of emb intolerance
  tm[lit > 6 & lit < 13 & cs == "TM2", "TMV+TM2"] <- 
    r$p_vis * 
    (1 - p_death2[cs == "TM2" & lit > 6 & lit < 13])
  ### Probability of death
  tm[lit > 6 & lit < 13 & cs == "TM2", "Death"] <- 
    p_death2[cs == "TM2" & lit > 6 & lit < 13] * 
    (1 - r$p_vis)
  ### Probability of remaining on treatment
  tm[lit > 6 & lit < 13 & cs == "TM2", "TM2"] <-  
    (1 - p_death2[cs == "TM2" & lit > 6 & lit < 13]) * 
    (1 - r$p_vis)
  
  ## TMV
  ### Probability of death
  tm[lit > 6 & lit < 13 & cs == "TMV", "Death"] <- 
    p_death2[cs == "TMV" & lit > 6 & lit < 13]
  ### Probability of remaining on treatment
  tm[lit > 6 & lit < 13 & cs == "TMV", "TMV"] <-
    (1 - p_death2[cs == "TMV" & lit > 6 & lit < 13])
  
  ## TM1+2
  ### Probability of emb intolerance
  tm[lit > 6 & lit < 13 & cs == "TM1+2", "TMV+TM1+2"] <- 
    r$p_vis * 
    (1 - p_death2[cs == "TM1+2" & lit > 6 & lit < 13])
  ### Probability of death
  tm[lit > 6 & lit < 13 & cs == "TM1+2", "Death"] <- 
    p_death2[cs == "TM1+2" & lit > 6 & lit < 13] * 
    (1 - r$p_vis)
  ### Probability of remaining on treatment
  tm[lit > 6 & lit < 13 & cs == "TM1+2", "TM1+2"] <- 
    (1 - p_death2[cs == "TM1+2" & lit > 6 & lit < 13]) * 
    (1 - r$p_vis)
  
  ## TMV+TM1
  ### Probability of death
  tm[lit > 6 & lit < 13 & cs == "TMV+TM1", "Death"] <- 
    p_death2[cs == "TMV+TM1" & lit > 6 & lit < 13]
  ### Probability of remaining on treatment
  tm[lit > 6 & lit < 13 & cs == "TMV+TM1", "TMV+TM1"] <- 
    (1 - p_death2[cs == "TMV+TM1" & lit > 6 & lit < 13])
  
  ## TMV+TM2
  ### Probability of death
  tm[lit > 6 & lit < 13 & cs == "TMV+TM2", "Death"] <- 
    p_death2[cs == "TMV+TM2" & lit > 6 & lit < 13]
  ### Probability of remaining on treatment
  tm[lit > 6 & lit < 13 & cs == "TMV+TM2", "TMV+TM2"] <- 
    (1 - p_death2[cs == "TMV+TM2" & lit > 6 & lit < 13])
  
  ## TMV+TM1+2
  ### Probability of death
  tm[lit > 6 & lit < 13 & cs == "TMV+TM1+2", "Death"] <- 
    p_death2[cs == "TMV+TM1+2" & lit > 6 & lit < 13]
  ### Probability of remaining on treatment
  tm[lit > 6 & lit < 13 & cs == "TMV+TM1+2", "TMV+TM1+2"] <- 
    (1 - p_death2[cs == "TMV+TM1+2" & lit > 6 & lit < 13])
  
  
  ##############################################################################
  # Month 13-18
  ##############################################################################
  ## On standard treatment
  ### Probability of emb intolerance
  tm[lit > 12 & lit < 19 & cs == "T", "TMV"] <- 
    r$p_vis * 
    (1 - r$p_cure_6mo) * 
    (1 - p_death2[cs == "T" & lit > 12 & lit < 19])
  ### Probability of cure
  tm[lit > 12 & lit < 19 & cs == "T", "Cure"] <- 
    r$p_cure_6mo * 
    (1 - r$p_vis) * 
    (1 - p_death2[cs == "T" & lit > 12 & lit < 19])
  ### Probability of death
  tm[lit > 12 & lit < 19 & cs == "T", "Death"] <- 
    p_death2[cs == "T" & lit > 12 & lit < 19] * 
    (1 - r$p_cure_6mo) * 
    (1 - r$p_vis)
  ### Probability of remaining on treatment
  tm[lit > 12 & lit < 19 & cs == "T", "T"] <-  
    (1 - p_death2[cs == "T" & lit > 12 & lit < 19]) * 
    (1 - r$p_cure_6mo) * 
    (1 - r$p_vis)
  
  ## On TM1
  ### Probability of emb intolerance
  tm[lit > 12 & lit < 19 & cs == "TM1", "TMV+TM1"] <- 
    r$p_vis * 
    (1 - r$p_cure_wo_rif_6mo) * 
    (1 - p_death2[cs == "TM1" & lit > 12 & lit < 19])
  ### Probability of cure
  tm[lit > 12 & lit < 19 & cs == "TM1", "Cure"] <- 
    r$p_cure_wo_rif_6mo * 
    (1 - r$p_vis) * 
    (1 - p_death2[cs == "TM1" & lit > 12 & lit < 19])
  ### Probability of death
  tm[lit > 12 & lit < 19 & cs == "TM1", "Death"] <- 
    p_death2[cs == "TM1" & lit > 12 & lit < 19] * 
    (1 - r$p_cure_wo_rif_6mo) * 
    (1 - r$p_vis)
  ### Probability of remaining on treatment
  tm[lit > 12 & lit < 19 & cs == "TM1", "TM1"] <- 
    (1 - p_death2[cs == "TM1" & lit > 12 & lit < 19]) * 
    (1 - r$p_cure_wo_rif_6mo) * 
    (1 - r$p_vis)
  
  ## On TM2
  ### Probability of emb intolerance
  tm[lit > 12 & lit < 19 & cs == "TM2", "TMV+TM2"] <- 
    r$p_vis * 
    (1 - r$p_cure_wo_azi_6mo) * 
    (1 - p_death2[cs == "TM2" & lit > 12 & lit < 19])
  ### Probability of cure
  tm[lit > 12 & lit < 19 & cs == "TM2", "Cure"] <- 
    r$p_cure_wo_azi_6mo * 
    (1 - r$p_vis) * 
    (1 - p_death2[cs == "TM2" & lit > 12 & lit < 19])
  ### Probability of death
  tm[lit > 12 & lit < 19 & cs == "TM2", "Death"] <- 
    p_death2[cs == "TM2" & lit > 12 & lit < 19] * 
    (1 - r$p_cure_wo_azi_6mo) * 
    (1 - r$p_vis)
  ### Probability of remaining on treatment
  tm[lit > 12 & lit < 19 & cs == "TM2", "TM2"] <- 
    (1 - p_death2[cs == "TM2" & lit > 12 & lit < 19]) * 
    (1 - r$p_cure_wo_azi_6mo) * 
    (1 - r$p_vis)
  
  ## On TMV
  ### Probability of cure
  tm[lit > 12 & lit < 19 & cs == "TMV", "Cure"] <- 
    r$p_cure_wo_emb_6mo * 
    (1 - p_death2[cs == "TMV" & lit > 12 & lit < 19])
  ### Probability of death
  tm[lit > 12 & lit < 19 & cs == "TMV", "Death"] <- 
    p_death2[cs == "TMV" & lit > 12 & lit < 19] * 
    (1 - r$p_cure_wo_emb_6mo)
  ### Probability of remaining on treatment
  tm[lit > 12 & lit < 19 & cs == "TMV", "TMV"] <-
    (1 - p_death2[cs == "TMV" & lit > 12 & lit < 19]) * 
    (1 - r$p_cure_wo_emb_6mo)
  
  ## On TM1+2
  ### Probability of emb intolerance
  tm[lit > 12 & lit < 19 & cs == "TM1+2", "TMV+TM1+2"] <- 
    r$p_vis * 
    (1 - min(r$p_cure_wo_rif_6mo, r$p_cure_wo_azi_6mo)) * 
    (1 - p_death2[cs == "TM1+2" & lit > 12 & lit < 19])
  ### Probability of cure
  tm[lit > 12 & lit < 19 & cs == "TM1+2", "Cure"] <- 
    min(r$p_cure_wo_rif_6mo, r$p_cure_wo_azi_6mo) * 
    (1 - p_death2[cs == "TM1+2" & lit > 12 & lit < 19]) * 
    (1 - r$p_vis)
  ### Probability of death
  tm[lit > 12 & lit < 19 & cs == "TM1+2", "Death"] <- 
    p_death2[cs == "TM1+2" & lit > 12 & lit < 19] * 
    (1 - min(r$p_cure_wo_rif_6mo, r$p_cure_wo_azi_6mo)) * 
    (1 - r$p_vis)
  ### Probability of remaining on treatment
  tm[lit > 12 & lit < 19 & cs == "TM1+2", "TM1+2"] <-
    (1 - p_death2[cs == "TM1+2" & lit > 12 & lit < 19]) * 
    (1 - min(r$p_cure_wo_rif_6mo, r$p_cure_wo_azi_6mo)) * 
    (1 - r$p_vis)
  
  ## On TMV+TM1
  ### Probability of cure
  tm[lit > 12 & lit < 19 & cs == "TMV+TM1", "Cure"] <- 
    min(r$p_cure_wo_rif_6mo, r$p_cure_wo_emb_6mo) * 
    (1 - p_death2[cs == "TMV+TM1" & lit > 12 & lit < 19])
  ### Probability of death
  tm[lit > 12 & lit < 19 & cs == "TMV+TM1", "Death"] <- 
    p_death2[cs == "TMV+TM1" & lit > 12 & lit < 19] * 
    (1 - min(r$p_cure_wo_rif_6mo, r$p_cure_wo_emb_6mo))
  ### Probability of remaining on treatment
  tm[lit > 12 & lit < 19 & cs == "TMV+TM1", "TMV+TM1"] <-
    (1 - p_death2[cs == "TMV+TM1" & lit > 12 & lit < 19]) * 
    (1 - min(r$p_cure_wo_rif_6mo, r$p_cure_wo_emb_6mo))
  
  ## On TMV+TM2
  ### Probability of cure
  tm[lit > 12 & lit < 19 & cs == "TMV+TM2", "Cure"] <- 
    min(r$p_cure_wo_emb_6mo, r$p_cure_wo_azi_6mo) * 
    (1 - p_death2[cs == "TMV+TM2" & lit > 12 & lit < 19])
  ### Probability of death
  tm[lit > 12 & lit < 19 & cs == "TMV+TM2", "Death"] <- 
    p_death2[cs == "TMV+TM2" & lit > 12 & lit < 19] * 
    (1 - min(r$p_cure_wo_emb_6mo, r$p_cure_wo_azi_6mo))
  ### Probability of remaining on treatment
  tm[lit > 12 & lit < 19 & cs == "TMV+TM2", "TMV+TM2"] <-
    (1 - p_death2[cs == "TMV+TM2" & lit > 12 & lit < 19]) * 
    (1 - min(r$p_cure_wo_emb_6mo, r$p_cure_wo_azi_6mo))
  
  ## On TMV+TM1+2
  ### Probability of cure
  tm[lit > 12 & lit < 19 & cs == "TMV+TM1+2", "Cure"] <- 
    min(r$p_cure_wo_emb_6mo, r$p_cure_wo_azi_6mo, r$p_cure_wo_rif_6mo) * 
    (1 - p_death2[cs == "TMV+TM1+2" & lit > 12 & lit < 19])
  ### Probability of death
  tm[lit > 12 & lit < 19 & cs == "TMV+TM1+2", "Death"] <- 
    p_death2[cs == "TMV+TM1+2" & lit > 12 & lit < 19] * 
    (1 - min(r$p_cure_wo_emb_6mo, r$p_cure_wo_azi_6mo, r$p_cure_wo_rif_6mo))
  ### Probability of remaining on treatment
  tm[lit > 12 & lit < 19 & cs == "TMV+TM1+2", "TMV+TM1+2"] <- 
    (1 - p_death2[cs == "TMV+TM1+2" & lit > 12 & lit < 19]) * 
    (1 - min(r$p_cure_wo_emb_6mo, r$p_cure_wo_azi_6mo, r$p_cure_wo_rif_6mo))
  ##############################################################################
  
  ##############################################################################
  # Month 19-23
  ##############################################################################
  ## On standard treatment
  ### Probability of emb intolerance
  tm[lit > 18 & lit < 24 & cs == "T", "TMV"] <- 
    r$p_vis * 
    (1 - r$p_cure_late) * 
    (1 - p_death2[cs == "T" & lit > 18 & lit < 24])
  ### Probability of cure
  tm[lit > 18 & lit < 24 & cs == "T", "Cure"] <- 
    r$p_cure_late * 
    (1 - r$p_vis) * 
    (1 - p_death2[cs == "T" & lit > 18 & lit < 24])
  ### Probability of death
  tm[lit > 18 & lit < 24 & cs == "T", "Death"] <- 
    p_death2[cs == "T" & lit > 18 & lit < 24] * 
    (1 - r$p_cure_late) * 
    (1 - r$p_vis)
  ### Probability of remaining on treatment
  tm[lit > 18 & lit < 24 & cs == "T", "T"] <-
    (1 - p_death2[cs == "T" & lit > 18 & lit < 24]) * 
    (1 - r$p_cure_late) * 
    (1 - r$p_vis)
  
  ## On TM1
  ### Probability of emb intolerance
  tm[lit > 18 & lit < 24 & cs == "TM1", "TMV+TM1"] <- 
    r$p_vis * 
    (1 - r$p_cure_wo_rif_late) * 
    (1 - p_death2[cs == "TM1" & lit > 18 & lit < 24])
  ### Probability of cure
  tm[lit > 18 & lit < 24 & cs == "TM1", "Cure"] <- 
    r$p_cure_wo_rif_late * 
    (1 - r$p_vis) * 
    (1 - p_death2[cs == "TM1" & lit > 18 & lit < 24])
  ### Probability of death
  tm[lit > 18 & lit < 24 & cs == "TM1", "Death"] <- 
    p_death2[cs == "TM1" & lit > 18 & lit < 24] * 
    (1 - r$p_cure_wo_rif_late) * 
    (1 - r$p_vis)
  ### Probability of remaining on treatment
  tm[lit > 18 & lit < 24 & cs == "TM1", "TM1"] <- 
    (1 - p_death2[cs == "TM1" & lit > 18 & lit < 24]) * 
    (1 - r$p_cure_wo_rif_late) * 
    (1 - r$p_vis)
  
  ## On TM2
  ### Probability of emb intolerance
  tm[lit > 18 & lit < 24 & cs == "TM2", "TMV+TM2"] <- 
    r$p_vis * 
    (1 - r$p_cure_wo_azi_late) * 
    (1 - p_death2[cs == "TM2" & lit > 18 & lit < 24])
  ### Probability of cure
  tm[lit > 18 & lit < 24 & cs == "TM2", "Cure"] <- 
    r$p_cure_wo_azi_late * 
    (1 - r$p_vis) * 
    (1 - p_death2[cs == "TM2" & lit > 18 & lit < 24])
  ### Probability of death
  tm[lit > 18 & lit < 24 & cs == "TM2", "Death"] <- 
    p_death2[cs == "TM2" & lit > 18 & lit < 24] * 
    (1 - r$p_cure_wo_azi_late) * 
    (1 - r$p_vis)
  ### Probability of remaining on treatment
  tm[lit > 18 & lit < 24 & cs == "TM2", "TM2"] <- 
    (1 - p_death2[cs == "TM2" & lit > 18 & lit < 24]) * 
    (1 - r$p_cure_wo_azi_late) * 
    (1 - r$p_vis)
  
  ## On TMV
  ### Probability of cure
  tm[lit > 18 & lit < 24 & cs == "TMV", "Cure"] <- 
    r$p_cure_wo_emb_late * 
    (1 - p_death2[cs == "TMV" & lit > 18 & lit < 24])
  ### Probability of death
  tm[lit > 18 & lit < 24 & cs == "TMV", "Death"] <- 
    p_death2[cs == "TMV" & lit > 18 & lit < 24] * 
    (1 - r$p_cure_wo_emb_late)
  ### Probability of remaining on treatment
  tm[lit > 18 & lit < 24 & cs == "TMV", "TMV"] <-
    (1 - p_death2[cs == "TMV" & lit > 18 & lit < 24]) * 
    (1 - r$p_cure_wo_emb_late)
  
  ## On TM1+2
  ### Probability of emb intolerance
  tm[lit > 18 & lit < 24 & cs == "TM1+2", "TMV+TM1+2"] <- 
    r$p_vis * 
    (1 - min(r$p_cure_wo_rif_late, r$p_cure_wo_azi_late)) * 
    (1 - p_death2[cs == "TM1+2" & lit > 18 & lit < 24])
  ### Probability of cure
  tm[lit > 18 & lit < 24 & cs == "TM1+2", "Cure"] <- 
    min(r$p_cure_wo_rif_late, r$p_cure_wo_azi_late) * 
    (1 - r$p_vis) * 
    (1 - p_death2[cs == "TM1+2" & lit > 18 & lit < 24])
  ### Probability of death
  tm[lit > 18 & lit < 24 & cs == "TM1+2", "Death"] <- 
    p_death2[cs == "TM1+2" & lit > 18 & lit < 24] * 
    (1 - min(r$p_cure_wo_rif_late, r$p_cure_wo_azi_late)) * 
    (1 - r$p_vis)
  ### Probability of remaining on treatment
  tm[lit > 18 & lit < 24 & cs == "TM1+2", "TM1+2"] <-
    (1 - p_death2[cs == "TM1+2" & lit > 18 & lit < 24]) * 
    (1 - min(r$p_cure_wo_rif_late, r$p_cure_wo_azi_late)) * 
    (1 - r$p_vis)
  
  ## On TMV+TM1
  ### Probability of cure
  tm[lit > 18 & lit < 24 & cs == "TMV+TM1", "Cure"] <- 
    min(r$p_cure_wo_rif_late, r$p_cure_wo_emb_late) * 
    (1 - p_death2[cs == "TMV+TM1" & lit > 18 & lit < 24])
  ### Probability of death
  tm[lit > 18 & lit < 24 & cs == "TMV+TM1", "Death"] <- 
    p_death2[cs == "TMV+TM1" & lit > 18 & lit < 24] * 
    (1 - min(r$p_cure_wo_rif_late, r$p_cure_wo_emb_late))
  ### Probability of remaining on treatment
  tm[lit > 18 & lit < 24 & cs == "TMV+TM1", "TMV+TM1"] <-
    (1 - p_death2[cs == "TMV+TM1" & lit > 18 & lit < 24]) * 
    (1 - min(r$p_cure_wo_rif_late, r$p_cure_wo_emb_late))
  
  ## On TMV+TM2
  ### Probability of cure
  tm[lit > 18 & lit < 24 & cs == "TMV+TM2", "Cure"] <- 
    min(r$p_cure_wo_emb_late, r$p_cure_wo_azi_late) * 
    (1 - p_death2[cs == "TMV+TM2" & lit > 18 & lit < 24])
  ### Probability of death
  tm[lit > 18 & lit < 24 & cs == "TMV+TM2", "Death"] <- 
    p_death2[cs == "TMV+TM2" & lit > 18 & lit < 24] * 
    (1 - min(r$p_cure_wo_emb_late, r$p_cure_wo_azi_late))
  ### Probability of remaining on treatment
  tm[lit > 18 & lit < 24 & cs == "TMV+TM2", "TMV+TM2"] <- 
    (1 - p_death2[cs == "TMV+TM2" & lit > 18 & lit < 24]) * 
    (1 - min(r$p_cure_wo_emb_late, r$p_cure_wo_azi_late))
  
  ## On TMV+TM1+2
  ### Probability of cure
  tm[lit > 18 & lit < 24 & cs == "TMV+TM1+2", "Cure"] <- 
    min(r$p_cure_wo_emb_late, r$p_cure_wo_azi_late,r$p_cure_wo_rif_late) * 
    (1 - p_death2[cs == "TMV+TM1+2" & lit > 18 & lit < 24])
  ### Probability of death
  tm[lit > 18 & lit < 24 & cs == "TMV+TM1+2", "Death"] <- 
    p_death2[cs == "TMV+TM1+2" & lit > 18 & lit < 24] * 
    (1 - min(r$p_cure_wo_emb_late, r$p_cure_wo_azi_late, r$p_cure_wo_rif_late))
  ### Probability of remaining on treatment
  tm[lit > 18 & lit < 24 & cs == "TMV+TM1+2", "TMV+TM1+2"] <-
    (1 - p_death2[cs == "TMV+TM1+2" & lit > 18 & lit < 24]) * 
    (1 - min(r$p_cure_wo_emb_late, r$p_cure_wo_azi_late, r$p_cure_wo_rif_late))
  ##############################################################################
  
  ##############################################################################
  # Month 24
  ##############################################################################
  ## On standard treatment
  ### Probability of cure
  tm[lit == 24 & cs == "T", "Cure"] <- 
    r$p_cure_late * 
    (1 - p_death2[cs == "T" & lit == 24])
  ### Probability of death
  tm[lit == 24 & cs == "T", "Death"] <- 
    p_death2[cs == "T" & lit == 24] * 
    (1 - r$p_cure_late)
  ### Probability of going to observation uncured
  tm[lit == 24 & cs == "T", "Obs"] <- 
    (1 - p_death2[cs == "T" & lit == 24]) * 
    (1 - r$p_cure_late)
  
  ## On TM1
  ### Probability of cure
  tm[lit == 24 & cs == "TM1", "Cure"] <- 
    r$p_cure_wo_rif_late * 
    (1 - p_death2[cs == "TM1" & lit == 24])
  ### Probability of death
  tm[lit == 24 & cs == "TM1", "Death"] <- 
    p_death2[cs == "TM1" & lit == 24] * 
    (1 - r$p_cure_wo_rif_late)
  ### Probability of going to observation uncured
  tm[lit == 24 & cs == "TM1", "Obs"] <- 
    (1 - p_death2[cs == "TM1" & lit == 24]) * 
    (1 - r$p_cure_wo_rif_late)
  
  ## On TM2
  ### Probability of cure
  tm[lit == 24 & cs == "TM2", "Cure"] <- 
    r$p_cure_wo_azi_late * 
    (1 - p_death2[cs == "TM2" & lit == 24])
  ### Probability of death
  tm[lit == 24 & cs == "TM2", "Death"] <- 
    p_death2[cs == "TM2" & lit == 24] * 
    (1 - r$p_cure_wo_azi_late)
  ### Probability of going to observation uncured
  tm[lit == 24 & cs == "TM2", "Obs"] <- 
    (1 - p_death2[cs == "TM2" & lit == 24]) * 
    (1 - r$p_cure_wo_azi_late)
  
  ## On TMV
  ### Probability of cure
  tm[lit == 24 & cs == "TMV", "Cure"] <- 
    r$p_cure_wo_emb_late * 
    (1 - p_death2[cs == "TMV" & lit == 24])
  ### Probability of death
  tm[lit == 24 & cs == "TMV", "Death"] <- 
    p_death2[cs == "TMV" & lit == 24] * 
    (1 - r$p_cure_wo_emb_late)
  ### Probability of going to observation uncured
  tm[lit == 24 & cs == "TMV", "Obs"] <- 
    (1 - p_death2[cs == "TMV" & lit == 24]) * 
    (1 - r$p_cure_wo_emb_late)
  
  ## On TM1+2
  ### Probability of cure
  tm[lit == 24 & cs == "TM1+2", "Cure"] <- 
    min(r$p_cure_wo_rif_late, r$p_cure_wo_azi_late) * 
    (1 - p_death2[cs == "TM1+2" & lit == 24])
  ### Probability of death
  tm[lit == 24 & cs == "TM1+2", "Death"] <- 
    p_death2[cs == "TM1+2" & lit == 24] * 
    (1 - min(r$p_cure_wo_rif_late, r$p_cure_wo_azi_late))
  ### Probability of going to observation uncured
  tm[lit == 24 & cs == "TM1+2", "Obs"] <- 
    (1 - p_death2[cs == "TM1+2" & lit == 24]) * 
    (1 - min(r$p_cure_wo_rif_late, r$p_cure_wo_azi_late))
  
  ## On TMV+TM1
  ### Probability of cure
  tm[lit == 24 & cs == "TMV+TM1", "Cure"] <- 
    min(r$p_cure_wo_rif_late, r$p_cure_wo_emb_late) * 
    (1 - p_death2[cs == "TMV+TM1" & lit == 24])
  ### Probability of death
  tm[lit == 24 & cs == "TMV+TM1", "Death"] <- 
    p_death2[cs == "TMV+TM1" & lit == 24] * 
    (1 - min(r$p_cure_wo_rif_late, r$p_cure_wo_emb_late))
  ### Probability of going to observation uncured
  tm[lit == 24 & cs == "TMV+TM1", "Obs"] <-
    (1 - p_death2[cs == "TMV+TM1" & lit == 24]) * 
    (1 - min(r$p_cure_wo_rif_late, r$p_cure_wo_emb_late))
  
  ## On TMV+TM2
  ### Probability of cure
  tm[lit == 24 & cs == "TMV+TM2", "Cure"] <- 
    min(r$p_cure_wo_emb_late, r$p_cure_wo_azi_late) * 
    (1 - p_death2[cs == "TMV+TM2" & lit == 24])
  ### Probability of going to observation uncured
  tm[lit == 24 & cs == "TMV+TM2", "Death"] <- 
    p_death2[cs == "TMV+TM2" & lit == 24] * 
    (1 - min(r$p_cure_wo_emb_late, r$p_cure_wo_azi_late))
  ### Probability of remaining not cured
  tm[lit == 24 & cs == "TMV+TM2", "Obs"] <- 
    (1 - p_death2[cs == "TMV+TM2" & lit == 24]) * 
    (1 - min(r$p_cure_wo_emb_late, r$p_cure_wo_azi_late))
  
  ## On TMV+TM1+2
  ### Probability of cure
  tm[lit == 24 & cs == "TMV+TM1+2", "Cure"] <- 
    min(r$p_cure_wo_emb_late, r$p_cure_wo_azi_late, r$p_cure_wo_rif_late) * 
    (1 - p_death2[cs == "TMV+TM1+2" & lit == 24])
  ### Probability of death
  tm[lit == 24 & cs == "TMV+TM1+2", "Death"] <- 
    p_death2[cs == "TMV+TM1+2" & lit == 24] * 
    (1 - min(r$p_cure_wo_emb_late, r$p_cure_wo_azi_late, r$p_cure_wo_rif_late))
  ### Probability of going to observation uncured
  tm[lit == 24 & cs == "TMV+TM1+2", "Obs"] <-
    (1 - p_death2[cs == "TMV+TM1+2" & lit == 24]) * 
    (1 - min(r$p_cure_wo_emb_late, r$p_cure_wo_azi_late, r$p_cure_wo_rif_late))
  
  if(any(lit == 25)) {
    stop(paste("Error in treatment lengths for individuals: "), 
          which(lit == 25))
  }
  ##############################################################################
  
  ## Observation
  tm[cs == "Obs" & !prev_vis_se, "T"] <- 
    r$p_prog * 
    (1 - p_death2[cs == "Obs" & !prev_vis_se]) * 
    (1 - r$p_spont_cure)
  tm[cs == "Obs" & prev_vis_se, "TMV"] <- 
    r$p_prog * 
    (1 - p_death2[cs == "Obs" & prev_vis_se]) * 
    (1 - r$p_spont_cure)
  tm[cs == "Obs", "Cure"] <- 
    r$p_spont_cure * 
    (1 - p_death2[cs == "Obs"]) * 
    (1 - r$p_prog)
  tm[cs == "Obs", "Death"] <- 
    p_death2[cs == "Obs"] * 
    (1 - r$p_spont_cure) * 
    (1 - r$p_prog)
  tm[cs == "Obs", "Obs"] <- 
    (1 - p_death2[cs == "Obs"]) * 
    (1 - r$p_spont_cure) * 
    (1 - r$p_prog)
  
  ## Cure
  tm[cs == "Cure", "Obs"] <- 
    r$p_recurrence * 
    (1 - p_death2[cs == "Cure"])
  tm[cs == "Cure", "Death"] <- 
    p_death2[cs == "Cure"] * 
    (1 - r$p_recurrence)
  tm[cs == "Cure", "Cure"] <- 
    (1 - p_death2[cs == "Cure"]) * 
    (1 - r$p_recurrence)
  
  return(tm)
}