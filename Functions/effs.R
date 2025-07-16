effs <- function(state_matrix, cy, r, prev_cc6m, 
                 prev_vis_se, obsTcounter, lit, util) {
  
  ## r = l_params
  ## Pass effs function to have input of vector of states for all individuals
  ## current state for all individuals + another cycle for all states in the cyclee before
  ## logical statememnt in vectorized form
  
  # list2env(l_params, environment(effs))
  cur_state <- state_matrix[, cy + 1]
  prev_state <- state_matrix[, cy]
  
  # if(cy == 600) {
  #   browser()
  # }
  
  util[cur_state == "T" | 
         cur_state == "TMV" |
         cur_state == "TM1" | 
         cur_state == "TM2" | 
         cur_state == "TM1+2" | 
         cur_state == "TMV+TM1" |
         cur_state == "TMV+TM2" | 
         cur_state == "TMV+TM1+2"] <- r$u_T  
  util[(prev_state == "T" & cur_state == "TMV") |
         (prev_state == "TM1" & cur_state == "TMV+TM1") |
         (prev_state == "TM2" & cur_state == "TMV+TM2") |
         (prev_state == "TM1+2" & cur_state == "TMV+TM1+2")] <- 
    r$U_SE_Etb - r$DisU_T
  util[(prev_state == "T" & cur_state == "TM1") |
         (prev_state == "TMV" & cur_state == "TMV+TM1") | 
         (prev_state == "TMV+TM2" & cur_state == "TMV+TM1+2")] <- 
    r$u_T - r$DisU_SE_Rif - r$DisU_T
  util[(prev_state == "T" & cur_state == "TM2") |
         (prev_state == "TMV" & cur_state == "TMV+TM2") |
         (prev_state == "TMV+TM1" & cur_state == "TMV+TM1+2")] <- 
    r$u_T - r$DisU_SE_Azi - r$DisU_T
  
  util[cur_state == "Obs"] <- r$u_obv
  util[cur_state == "Cure"] <- r$u_cure
  util[cur_state == "Death"] <- 0
  
  if(any(is.na(util))) { stop("Error in utility logic") }
  
  # Updating counter for progression if required treatment
  obsTcounter[((prev_state == "Obs") & (cur_state == "T")) | 
                ((prev_state == "Obs") & (cur_state == "TMV"))] <- 
    obsTcounter[((prev_state == "Obs") & (cur_state == "T")) |
                  ((prev_state == "Obs") & (cur_state == "TMV"))] + 1
  # Updating counter for cc6m
  prev_cc6m[lit <= 18 & cur_state == "Cure" &
              (prev_state == "T" | 
                 prev_state == "TMV" |
                 prev_state == "TM1" | 
                 prev_state == "TM2" | 
                 prev_state == "TM1+2" | 
                 prev_state == "TMV+TM1" |
                 prev_state == "TMV+TM2" | 
                 prev_state == "TMV+TM1+2")] <- T
  # Updating counter for prev_vis_se
  prev_vis_se[cur_state == "TMV" | 
                cur_state == "TMV+TM1" | 
                cur_state == "TMV+TM2" | 
                cur_state == "TMV+TM1+2"] <- T
  # Updating length in treatment
  lit[(cur_state == "T" | 
         cur_state == "TMV" |
         cur_state == "TM1" | 
         cur_state == "TM2" | 
         cur_state == "TM1+2" | 
         cur_state == "TMV+TM1" |
         cur_state == "TMV+TM2" | 
         cur_state == "TMV+TM1+2") &
        (prev_state == "T" | 
           prev_state == "TMV" |
           prev_state == "TM1" | 
           prev_state == "TM2" | 
           prev_state == "TM1+2" | 
           prev_state == "TMV+TM1" |
           prev_state == "TMV+TM2" | 
           prev_state == "TMV+TM1+2")] <- lit[(cur_state == "T" | 
                                                 cur_state == "TMV" |
                                                 cur_state == "TM1" | 
                                                 cur_state == "TM2" | 
                                                 cur_state == "TM1+2" | 
                                                 cur_state == "TMV+TM1" |
                                                 cur_state == "TMV+TM2" | 
                                                 cur_state == "TMV+TM1+2") &
                                                (prev_state == "T" | 
                                                   prev_state == "TMV" |
                                                   prev_state == "TM1" | 
                                                   prev_state == "TM2" | 
                                                   prev_state == "TM1+2" | 
                                                   prev_state == "TMV+TM1" |
                                                   prev_state == "TMV+TM2" | 
                                                   prev_state == "TMV+TM1+2")] + 1
  lit[!((cur_state == "T" | 
          cur_state == "TMV" |
          cur_state == "TM1" | 
          cur_state == "TM2" | 
          cur_state == "TM1+2" | 
          cur_state == "TMV+TM1" |
          cur_state == "TMV+TM2" | 
          cur_state == "TMV+TM1+2") &
         (prev_state == "T" | 
            prev_state == "TMV" |
            prev_state == "TM1" | 
            prev_state == "TM2" | 
            prev_state == "TM1+2" | 
            prev_state == "TMV+TM1" |
            prev_state == "TMV+TM2" | 
            prev_state == "TMV+TM1+2"))] <- 0
  
  # Disutility of progression if observed and required treatment
  util <- util - (obsTcounter * r$DisU_obsT)
  util[util < 0] <- 0
  
  return(
    list(
      "util" = util,
      "prev_cc6m" = prev_cc6m,
      "prev_vis_se" = prev_vis_se,
      "obsTcounter" = obsTcounter,
      "lit" = lit
    )
  )
  
}