library(tidyverse)
library(parallel)
source("microsim.R")
source("Functions/life_expectancy.R")
cycle_length <- 1 / 12

###########################################################
### Baseline values
ni <- 1e4   # Number of people
fl <- 40    # Years

# Obtaining baseline values
baseline.vals <- read.csv("values/baseline.csv", 
                          row.names = "variable_name")
params_base <- as.list(baseline.vals$value)
params_max <- as.list(baseline.vals$max_value)
params_min <- as.list(baseline.vals$min_value)

params_name <- rownames(baseline.vals)
names(params_base) <- params_name
names(params_min) <- params_name
names(params_max) <- params_name

# For DSA
n_dsa <- 10
dsa_val <- map2(.x = params_min, 
                .y = params_max, 
                \(x,y) seq(x, y, length.out = n_dsa))

# Reading life table from Canada Census for females in 2021
life.table <- read.csv("values/2022_lifetable.csv") |> 
  select("Age", "qx") |>  # qx is the death rate
  rename(age = Age, base_rate = qx) |> 
  # Changing age to be just the number
  mutate(age = as.numeric(str_split_i(age, " ", i = 1)))

# Retrieves the rate of death from the life table 
# (bcc = base culture conversion)
params_base$p_death <- life.table |> 
  # Adding additional rows so that we can align ages
  add_row(age = 111:120, base_rate = rep(1, 10)) |> 
  rename(control_rate = base_rate) |> 
  mutate(ntm_rate = control_rate * params_base$hr) |> 
  # Death rate for those who had culture conversion
  # in first 3 months
  mutate(cc6m_rate = ntm_rate * params_base$hr_cc6m) |> 
  as.matrix()

## For debugging. Comment out the next part in the working 
## model -------------------------------------------------
# list2env(params_base, .GlobalEnv)
# r <- params_base
# l_params <- params_base
# n_p = 100
# fup = 30
# cycle_length = 1 / 12
# given_seed = 100
# arm <- "treat"
# cycle <- 1

###########################################################
### Running the base model
set.seed(100)
b.t <- microsim_model(params_base, "treat", n_p = ni, fup = fl)
b.o <- microsim_model(params_base, "observe", n_p = ni, fup = fl)

base_result <- tibble(
  treatment = c("Treatment", "Observation"),
  avg = c(mean(b.t$util_total), mean(b.o$util_total)),
  se = c(sd(b.t$util_total), sd(b.o$util_total)) / sqrt(ni),
  life_expectancy = c(mean(life_expectancy(b.t$state_matrix)),
                      mean(life_expectancy(b.o$state_matrix)))) |> 
  mutate(across(2:4, ~ .x / 12))

# saveRDS(base_result, file = "Results/base.RData")
# base_result <- readRDS("Results/base.RData")

###########################################################
### Running the one way DSA
# not_in_dsa <- c("u_death")
not_in_dsa <- c("u_death", "discount")

# Initiate empty list
owsa_dsa <- lapply(dsa_val, \(x) 
                   tibble(value = x,
                          trt_util = NA, trt_le = NA,
                          obs_util = NA, obs_le = NA))
owsa_dsa <- owsa_dsa[!names(owsa_dsa) %in% not_in_dsa]
dsa_params <- names(owsa_dsa)

n_cores <- detectCores(logical = FALSE)
pb <- txtProgressBar(min = 1, max = length(dsa_params), style = 3)

for(i in 1:length(dsa_params)){
  # Seed
  s <- i + 2
  
  # Selecting parameter
  param <- dsa_params[i]
  
  # Preparing parameter sets for each variable to be varied
  params_tmp <- replicate(n_dsa, params_base)
  params_tmp[param, ] <- dsa_val[[param]]
  
  # Fix death if we are changing HR
  if(param == "hr") {
    pd <- params_tmp["p_death",][[1]] # all identical
    for(ii in 1:n_dsa){
      new_hr <- unlist(params_tmp["hr",])
      params_tmp["p_death",][[ii]][,"ntm_rate"] = 
        pd[,"control_rate"] * new_hr[ii]
      params_tmp["p_death",][[ii]][,"cc6m_rate"] = 
        params_tmp["p_death",][[ii]][,"ntm_rate"] * 
        params_base$hr_cc6m
    }
  } else if(param == "hr_cc6m") {
    pd <- params_tmp["p_death",][[1]] # all identical
    for(ii in 1:n_dsa){
      new_hrcc6m <- unlist(params_tmp["hr_cc6m",])
      params_tmp["p_death",][[ii]][,"cc6m_rate"] = 
        pd[,"ntm_rate"] * new_hrcc6m[[ii]]
    }
  }
  
  # Running Treatment
  set.seed(s)
  tmp.t <- mclapply(1:n_dsa, \(x) 
                    microsim_model(params_tmp[,x], arm = "treat",
                                   n_p = ni, fup = fl),
                    mc.cores = n_cores - 1)
  
  # Observation
  set.seed(s + 1)
  tmp.o <- mclapply(1:n_dsa, \(x) 
                    microsim_model(params_tmp[,x], arm = "observe",
                                   n_p = ni, fup = fl),
                    mc.cores = n_cores - 1)
  
  # Extracting values
  owsa_dsa[[param]]$trt_util <- 
    unlist(lapply(1:n_dsa, \(x) mean(tmp.t[[x]]$util_total)))
  owsa_dsa[[param]]$trt_le <- 
    unlist(lapply(1:n_dsa, \(x) 
                  mean(life_expectancy(tmp.t[[x]]$state_matrix))))
  
  owsa_dsa[[param]]$obs_util <- 
    unlist(lapply(1:n_dsa, \(x) mean(tmp.o[[x]]$util_total)))
  owsa_dsa[[param]]$obs_le <- 
    unlist(lapply(1:n_dsa, \(x) 
                  mean(life_expectancy(tmp.o[[x]]$state_matrix))))
  
  setTxtProgressBar(pb, i)
  
}


# saveRDS(owsa_dsa, file = "Results/owsa_dsa.RData")
# owsa_dsa <- readRDS("Results/owsa_dsa.RData")

# owsa_dsa |> 
#   bind_rows() |> 
#   mutate(util_trt_pref = trt_util > obs_util) |> 
#   mutate(util_le_pref = trt_le > obs_le) |> 
#   summarize(
#     util_trt_pref = all(util_trt_pref),
#     util_le_pref = all(util_le_pref)
#   )
##    # A tibble: 1 × 2
##    util_trt_pref util_le_pref
##    <lgl>         <lgl>       
##   1 TRUE          TRUE   

###########################################################
### Two way sensitivity analyses
# Initialize list
twsa_dsa <- list()

# Probability of spontaneous culture conversion and observation 
# to treatment
twsa_dsa[[1]] <- expand_grid(p_spont_cure = dsa_val$p_spont_cure, 
                             p_prog = dsa_val$p_prog)

# Utility of observation and utility of treatment
twsa_dsa[[2]] <- expand_grid(u_obv = dsa_val$u_obv, 
                             u_T = dsa_val$u_T)

# Utility of treatment and utility of cure
twsa_dsa[[3]] <- expand_grid(u_cure = dsa_val$u_cure, 
                             u_T = dsa_val$u_T)

# Utility of emb intolerance and ethambutol medication 
# intolerance
twsa_dsa[[4]] <- expand_grid(p_vis = dsa_val$p_vis, 
                             U_SE_Etb = dsa_val$U_SE_Etb)

# HR of culture conversion in 3 months of treatment and HR of 
# patients with MAC-PD
twsa_dsa[[5]] <- expand_grid(hr = dsa_val$hr, 
                             hr_cc6m = dsa_val$hr_cc6m)

twsa_dsa <- lapply(twsa_dsa, \(x) 
                   mutate(x, "trt_util" = NA, "trt_le" = NA, 
                          "obs_util" = NA, "obs_le" = NA))

# Sequence variables
move_by <- 10
s_seq <- seq(1, n_dsa ^ 2, by = move_by)

for(i in 1:length(twsa_dsa)) {
  cat("\nCombination", i, "of", length(twsa_dsa), "\n")
  seed <- (i * 2)
  
  # Names of variables included and excluded
  dsa_param_names <- colnames(twsa_dsa[[i]][1:2])
  other_params <- within(params_base, rm(list = dsa_param_names))
  
  # Deal with death separately
  other_params2 <- within(other_params, rm("p_death"))
  tmp_params <- lapply(other_params2, \(x) rep(x, move_by))
  tmp_params$p_death <- rep(list(other_params[["p_death"]]), move_by)
  
  for(ii in s_seq) {
    cat("\r  Part", ii, "-", ii + move_by - 1, 
        "of", dim(twsa_dsa[[i]])[1])
    
    # Make temporary df with values of sensitivity analysis
    tmp_seq <- ii:(ii + move_by - 1)
    tmp_df <- twsa_dsa[[i]][tmp_seq, 1:2]
    r <- tmp_params
    
    r[[dsa_param_names[1]]] <- tmp_df[[dsa_param_names[1]]]
    r[[dsa_param_names[2]]] <- tmp_df[[dsa_param_names[2]]]
    
    # Incorporating HR into death table
    if(any(dsa_param_names == "hr")) {
      for (iii in 1:length(r$p_death)) {
        r[["p_death"]][[iii]][,"ntm_rate"] <- 
          r[["p_death"]][[iii]][,"control_rate"] * r[["hr"]][iii]
        r[["p_death"]][[iii]][,"cc6m_rate"] <- 
          r[["p_death"]][[iii]][,"ntm_rate"] * r[["hr_cc6m"]][iii]
      }
    }
    
    # Change orientation
    r <- lapply(1:move_by, \(x) lapply(r, `[`, x))
    
    # Unlist death
    for(iii in 1:move_by) {
      r[[iii]][["p_death"]] <- r[[iii]][["p_death"]][[1]]
    }
    
    # Run the analysis, store in temporary `tmp` lists
    set.seed(seed)
    tmp.t <- mclapply(r, \(x) 
                      microsim_model(x, arm = "treat", n_p = ni, fup = fl),
                      mc.cores = n_cores - 1)
    
    tmp.o <- mclapply(r, \(x) 
                      microsim_model(x, arm = "observe", n_p = ni, fup = fl),
                      mc.cores = n_cores - 1)
    
    # Extract values from list into summary values
    twsa_dsa[[i]][["trt_util"]][tmp_seq] <-  
      unlist(lapply(1:move_by, \(x) mean(tmp.t[[x]]$util_total)))
    twsa_dsa[[i]][["trt_le"]][tmp_seq] <- 
      unlist(lapply(1:move_by, 
                    \(x) mean(life_expectancy(tmp.t[[x]]$state_matrix))))
    twsa_dsa[[i]][["obs_util"]][tmp_seq] <- 
      unlist(lapply(1:move_by, \(x) mean(tmp.o[[x]]$util_total)))
    twsa_dsa[[i]][["obs_le"]][tmp_seq] <- 
      unlist(lapply(1:move_by, 
                    \(x) mean(life_expectancy(tmp.o[[x]]$state_matrix))))
    
  }
}

# saveRDS(twsa_dsa, file = "Results/twsa_dsa.RData")
# twsa_dsa <- readRDS("Results/twsa_dsa.RData")

twsa_dsa_summary <- twsa_dsa |>
  map(~ mutate(.x,
               util_strategy = ifelse(trt_util > obs_util,
                                      "Treatment", "Observation"),
               le_strategy = ifelse(trt_le > obs_le,
                                    "Treatment", "Observation")))

all(bind_rows(twsa_dsa_summary)$util_strategy == "Treatment")
all(bind_rows(twsa_dsa_summary)$le_strategy == "Treatment")

## plot
twsa_plot1 <- twsa_dsa_summary[[5]] |>
  ggplot(aes(hr, hr_cc6m)) +
  geom_tile(aes(fill = le_strategy), colour = "black") +
  labs(x = "Hazard ratio of mortality for patients with MAC-PD",
       y = paste("Hazard ratio of martality for patients with\n",
                 "early culture conversion"),
       fill = "Strategy") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_brewer(palette = "Set2") +
  theme_linedraw()

twsa_plot2 <- twsa_dsa_summary[[5]] |>
  ggplot(aes(hr, hr_cc6m)) +
  geom_tile(aes(fill = util_strategy), colour = "black") +
  labs(x = "Hazard ratio of mortality for patients with MAC-PD",
       y = paste("Hazard ratio of martality for patients with\n",
                 "early culture conversion"),
       fill = "Strategy") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_brewer(palette = "Set1") +
  theme_linedraw()

###########################################################
### Plotting OWSA
owsa_dsa_long <- lapply(
  owsa_dsa, \(x)
  x |> 
    pivot_longer(
      cols = !value,
      names_to = c("strategy", "stat"),
      names_pattern = "(.*)_(.*)",
      values_to = "v"
    ) |> 
    pivot_wider(
      names_from = "stat",
      values_from = "v"
    )) |> 
  bind_rows(.id = "name") |> 
  mutate(
    strategy = case_when(
      strategy == "obs" ~ "Observation",
      strategy == "trt" ~ "Treatment"), 
    name = as_factor(name),
    # Change life expectancy to years
    le = le / 12) |> 
  # Add names for plotting
  left_join(rownames_to_column(baseline.vals, var = "name"), 
            by = "name") |>
  # Clean up
  select(!c("value.y", "min_value", "max_value")) |> 
  rename(value = value.x) |> 
  # Changing rate to yearly probability
  mutate(
    value = ifelse(str_starts(value, "p_"), 
                   r_p(p_r(value), 12),
                   value)
  ) |> 
  mutate(full_name = str_replace_all(full_name, "\\\\n", "\n")) |> 
  filter(name != "DisU_T")

owsa_plot1 <- owsa_dsa_long |> 
  ggplot(aes(value, util)) + 
  geom_line(aes(colour = strategy), linewidth = 1.5) +
  facet_wrap(vars(full_name), scales = "free") + 
  theme_linedraw() +
  labs(x = "Parameter Value", y = "QALYs", colour = "Strategy") + 
  scale_colour_brewer(palette = "Set1") +
  scale_y_continuous(
    labels = scales::number_format(accuracy = 0.1)) +
  theme(panel.grid = element_blank(),
        text = element_text(size = 14))

owsa_plot2 <- owsa_dsa_long |> 
  ggplot(aes(value, le)) + 
  geom_line(aes(colour = strategy), linewidth = 1.5) +
  facet_wrap(vars(full_name), scales = "free") + 
  theme_linedraw() +
  labs(x = "Parameter Value", y = "Life Exectancy (Years)", 
       colour = "Strategy") + 
  scale_colour_brewer(palette = "Set2") +
  scale_y_continuous(
    labels = scales::number_format(accuracy = 0.1)) +
  theme(panel.grid = element_blank(),
        text = element_text(size = 14))

owsa_plot2.1 <- owsa_dsa_long |>
  filter(name == "hr_cc6m") |> 
  ggplot(aes(value, le)) + 
  geom_line(aes(colour = strategy), linewidth = 1.5) +
  facet_wrap(vars(full_name), scales = "free") + 
  theme_linedraw() +
  labs(x = "Parameter Value", y = "Life Exectancy (Years)", 
       colour = "Strategy") + 
  scale_colour_brewer(palette = "Set2") +
  scale_y_continuous(
    labels = scales::number_format(accuracy = 0.1)) +
  theme(panel.grid = element_blank(),
        text = element_text(size = 14))

ggsave("owsa_plot1.pdf", plot = owsa_plot1, path = "Plots/",
       height = 12, width = 19, units = "in")
ggsave("owsa_plot2.pdf", plot = owsa_plot2, path = "Plots/",
       height = 12, width = 19, units = "in")
ggsave("owsa_plot2.1.pdf", plot = owsa_plot2.1, path = "Plots/",
       height = 4, width = 6, units = "in")
ggsave("twsa_plot1.pdf", plot = twsa_plot1, path = "Plots/",
       height = 4, width = 6, units = "in")



###########################################################
### Statistics from baseline probabilities (Table 2)

###
# For table
# Probabilities in yearly percent
baseline.vals |> 
  rownames_to_column(var = "name") |> 
  tibble() |> 
  filter(str_starts(name, "p_")) |> 
  filter(str_ends(name, "6mo") | str_ends(name, "late") | 
           str_starts(name, "p_trtd") | name == "p_vis") |> 
  mutate(across(2:4, ~ r_p(p_r(.x), 6) * 100))

baseline.vals |> 
  rownames_to_column(var = "name") |> 
  tibble() |> 
  filter(str_starts(name, "p_")) |> 
  filter(!(str_ends(name, "6mo") | str_ends(name, "late") | 
           str_starts(name, "p_trtd"))) |> 
  mutate(across(2:4, ~ r_p(p_r(.x), 12) * 100))
  
baseline.vals |> 
  rownames_to_column(var = "name") |> 
  tibble() |> 
  filter(!str_starts(name, "p_"))


# Ten-year survival
sum(life_expectancy(b.t$state_matrix[,1:(1 + 10 * 12)]) == (10 * 12)) / ni
# [1] 0.8773
sum(life_expectancy(b.o$state_matrix[,1:(1 + 10 * 12)]) == (10 * 12)) / ni
# [1] 0.8292

## Time on treatment
mean(rowSums(b.t$state_matrix == "T" |
               b.t$state_matrix == "TMV" | 
               b.t$state_matrix == "TM1" |
               b.t$state_matrix == "TM2" |
               b.t$state_matrix == "TM1+2" |
               b.t$state_matrix == "TMV+TM1" | 
               b.t$state_matrix == "TMV+TM2" |
               b.t$state_matrix == "TMV+TM1+2"))
# [1] 24.5726

mean(rowSums(b.o$state_matrix == "T" |
               b.o$state_matrix == "TMV" | 
               b.o$state_matrix == "TM1" |
               b.o$state_matrix == "TM2" |
               b.o$state_matrix == "TM1+2" |
               b.o$state_matrix == "TMV+TM1" | 
               b.o$state_matrix == "TMV+TM2" |
               b.o$state_matrix == "TMV+TM1+2"))
# [1] 7.8998

## Progress to treatment
sum(b.t$obsTcounter > 0) / ni
# [1] 0.3199
sum(b.o$obsTcounter > 0) / ni
# [1] 0.3612

## Number of individuals on observation requiring treatment
sum(rowSums(b.o$state_matrix == "T" |
              b.o$state_matrix == "TMV" | 
              b.o$state_matrix == "TM1" |
              b.o$state_matrix == "TM2" |
              b.o$state_matrix == "TM1+2" |
              b.o$state_matrix == "TMV+TM1" | 
              b.o$state_matrix == "TMV+TM2" |
              b.o$state_matrix == "TMV+TM1+2") > 0) / ni
# [1] 0.3612

## Individuals with rifampin intolerance
sum(rowSums(b.t$state_matrix == "TM1" | 
              b.t$state_matrix == "TM1+2" |
              b.t$state_matrix == "TMV+TM1" |
              b.t$state_matrix == "TMV+TM1+2") > 1) / ni
# [1] 0.0935
sum(rowSums(b.o$state_matrix == "TM1" | 
              b.o$state_matrix == "TM1+2" |
              b.o$state_matrix == "TMV+TM1" |
              b.o$state_matrix == "TMV+TM1+2") > 1) / ni
# [1] 0.0308

## Individuals with azi intolerance
sum(rowSums(b.t$state_matrix == "TM2" | 
              b.t$state_matrix == "TM1+2" |
              b.t$state_matrix == "TMV+TM2" |
              b.t$state_matrix == "TMV+TM1+2") > 1) / ni
# [1] 0.0406

sum(rowSums(b.o$state_matrix == "TM2" | 
              b.o$state_matrix == "TM1+2" |
              b.o$state_matrix == "TMV+TM2" |
              b.o$state_matrix == "TMV+TM1+2") > 1) / ni
# [1] 0.0121


## Individuals with ethambutol intolerance
sum(rowSums(b.t$state_matrix == "TMV" | 
              b.t$state_matrix == "TMV+TM1" |
              b.t$state_matrix == "TMV+TM2" |
              b.t$state_matrix == "TMV+TM1+2") > 1) / ni
# [1] 0.013
sum(rowSums(b.o$state_matrix == "TMV" | 
              b.o$state_matrix == "TMV+TM1" |
              b.o$state_matrix == "TMV+TM2" |
              b.o$state_matrix == "TMV+TM1+2") > 1) / ni
# [1] 0.0035

## More than one adverse event
sum(rowSums(b.t$state_matrix == "TM1+2" | 
              b.t$state_matrix == "TMV+TM1" |
              b.t$state_matrix == "TMV+TM2" |
              b.t$state_matrix == "TMV+TM1+2") > 1) / ni
# [1] 0.0041
sum(rowSums(b.o$state_matrix == "TM1+2" | 
              b.o$state_matrix == "TMV+TM1" |
              b.o$state_matrix == "TMV+TM2" |
              b.o$state_matrix == "TMV+TM1+2") > 1) / ni
# [1] 0.001

## Individuals in the cure state
sum(rowSums(b.t$state_matrix == "Cure") > 1) / ni
# [1] 0.9592
mean(rowSums(b.t$state_matrix == "Cure")) / 12
# [1] 11.60068
sum(rowSums(b.o$state_matrix == "Cure") > 1) / ni
# [1] 0.9016
mean(rowSums(b.o$state_matrix == "Cure")) / 12
# [1] 8.815425

