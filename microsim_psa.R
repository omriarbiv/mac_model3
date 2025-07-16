library(tidyverse)
library(parallel)
library(patchwork)
source("microsim.R")
source("Functions/gen_psa_params.R")
source("Functions/life_expectancy.R")


baseline.vals <- read.csv("values/baseline.csv", 
                          row.names = "variable_name")
params_base <- as.list(baseline.vals$value)
params_name <- rownames(baseline.vals)
names(params_base) <- params_name

# Adding information for hazard ratio
params_base$hrm1 <- log(baseline.vals["hr_cc6m", "value"])
params_base$hrv1 <- ((log(0.74) - log(baseline.vals["hr_cc6m", "min_value"])) 
                     / 3.92) ^ 2

params_base$hrm2 <- log(baseline.vals["hr", "value"])
params_base$hrv2 <- ((baseline.vals["hr", "max_value"] - 
                        baseline.vals["hr", "min_value"]) / 3.92) ^ 2

life.table <- read.csv("values/2022_lifetable.csv") |>
  select("Age", "qx") |>  # qx is the death rate
  rename(age = Age, control_rate = qx) |>
  # Changing age to be just the number
  mutate(age = as.numeric(str_split_i(age, " ", i = 1))) |>
  as.matrix()

n_cores <- detectCores(logical = FALSE)

run_psa <- function(iter, seed = NULL, save_psa = T,
                    params_base = params_base, 
                    baseline.vals = baseline.vals, 
                    life.table = life.table) {
  if(is.null(seed)) { seed <- 1000 }
  set.seed((seed + iter[length(iter)]))
  
  # Generating random parameters from their distributions
  d_psa_params <- gen_psa_params(n_psa = iter[length(iter)],
                                 seed = seed, pr = 0,
                                 psa_params = params_base,
                                 life.table = life.table) 
  # Re-organizing psa_params to be in the same style as params_base
  psa_params <- lapply(iter, \(x) lapply(d_psa_params, `[`, x))
  for(i in 1:length(iter)) {
    psa_params[[i]]$p_death <- matrix(
      unlist(psa_params[[i]]$p_death, use.names = T), 
      ncol = 4, byrow = F, 
      dimnames = list(NULL, c("age", "control_rate", "ntm_rate", "cc6m_rate"))
    )
  }

  # Running PSA for treatment & obs with random parameters
  print("Running Treatment Arms")
  psa.t <- mclapply(psa_params, \(x) 
               microsim_model(x, arm = "treat", n_p = ni, fup = fl, 
                              cycle_length = 1 / 12),
              mc.cores = n_cores - 1)
  print("Running Observational Arms")
  psa.o <- mclapply(psa_params, \(x) 
               microsim_model(x, arm = "observe",  n_p = ni, fup = fl, 
                              cycle_length = 1 / 12),
               mc.cores = n_cores - 1)
  
  # Combining parameters in matrix to save
  psa_result <- lapply(1:length(iter), \(x)
                       cbind(psa.t[[x]]$util_total, 
                             life_expectancy(psa.t[[x]]$state_matrix),
                             psa.o[[x]]$util_total,
                             life_expectancy(psa.o[[x]]$state_matrix)))
  
  for(i in 1:length(iter)) {
    colnames(psa_result[[i]]) <- c("trt_util", "trt_le",
                                   "obs_util", "obs_le")
  }
  
  if(save_psa) {
    saveRDS(psa_result, file = paste0("Results/psa_n", 
                                 iter[1], "_", iter[length(iter)], 
                                 ".RData"))
  }
  return(psa_result)
}

###########################################################
### Running PSA

ni <- 1e4
fl <- 40
niter <- 1e3
move_by <- 10
s_seq <- seq(1, niter, by = move_by)

for (i in 1:length(s_seq)) {
  start_time <- Sys.time()
  iter <- s_seq[i]:(s_seq[i] + move_by - 1)
  cat("Running sequence ", iter[1], " to ", 
               iter[length(iter)], "\n")
  pp2 <- run_psa(iter, params_base = params_base, 
                 life.table = life.table)
  print(Sys.time() - start_time)
}

###########################################################
### Analyzing PSA

# Find files in Results folder
result_files <- list.files("Results")
psa_files <- str_c("Results/", 
                   result_files[str_starts(result_files, "psa_")])

# Read in raw result and obtain means and deltas
psa_result <- map(psa_files, readRDS) |> 
  tibble() |> 
  rename("raw" = 1) |> 
  unnest(raw) |> 
  mutate(
    trt_util = map(raw, ~ mean(.x[, "trt_util"])),
    obs_util = map(raw, ~ mean(.x[, "obs_util"])),
    delta_util = map(raw, ~ mean(.x[, "trt_util"] - .x[, "obs_util"])),
    trt_le = map(raw, ~ mean(.x[, "trt_le"])),
    obs_le = map(raw, ~ mean(.x[, "obs_le"])),
    delta_le = map(raw, ~ mean(.x[, "trt_le"] - .x[, "obs_le"])),
  ) |> 
  unnest(!1) |> 
  mutate(across(!raw, ~ .x / 12))

# Summarize data
psa_summary <- psa_result |> 
  pivot_longer(
    cols = starts_with(c("trt", "obs", "delta")),
    names_to = c("strategy", ".value"),
    names_sep = "_"
  ) |> 
  mutate(
    avg_util = mean(util),
    avg_le = mean(le), 
    mcse_util = sd(util) / sqrt(niter),
    mcse_le = sd(le) / sqrt(niter),
    l95_util = quantile(util, 0.025),
    u95_util = quantile(util, 0.975),
    l95_le = quantile(le, 0.025),
    u95_le = quantile(le, 0.975),
    .by = strategy
  ) |> 
  select(!c(raw, util, le)) |>
  # relocate(starts_with("mcse"), .after = "delta_le") |> 
  slice(1:3) |>
  pivot_longer(
    cols = !1,
    names_to = c(".value", "outcome"),
    names_sep = "_"
  ) |> 
  arrange(desc(outcome), desc(strategy))

# Plot data
psa_plot_data <- psa_result |> 
  select(!raw) |> 
  pivot_longer(
    cols = ends_with(c("util", "le")),
    names_to = c("strategy", ".value"),
    names_sep = "_"
    ) |> 
  mutate(avg_util = mean(util), 
         avg_le = mean(le), 
         .by = "strategy") |> 
  mutate(
    strategy = case_when(strategy == "trt" ~ "Treatment",
                         strategy == "obs" ~ "Observation",
                         strategy == "delta" ~ "Treatment to Observation Delta"))

psa_plot1.1 <- psa_plot_data |> 
  filter(strategy != "Treatment to Observation Delta") |> 
  ggplot(aes(util)) + 
  geom_histogram(aes(fill = strategy, colour = strategy),
                 bins = 30, alpha = 0.8, position = "identity") +
  geom_vline(aes(xintercept = avg_util), linetype = 2) + 
  scale_fill_brewer(palette = "Set1") +
  scale_colour_brewer(palette = "Set1") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), limit = c(0, 200)) +
  facet_wrap(vars(strategy), nrow = 2) + 
  labs(x = "QALYs", y = "Count", fill = "Strategy") +
  theme_linedraw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

psa_plot1.2 <- psa_plot_data |>
  filter(strategy == "Treatment to Observation Delta") |> 
  # mutate(util0 = ifelse(util > 0, "G", "L")) |> 
  ggplot(aes(util)) + 
  geom_histogram(
    # aes(colour = util0, fill = util0),
    colour = "green4", fill = "green4",
    alpha = 0.8, bins = 30, 
    position = "identity") +
  # geom_vline(aes(xintercept = avg)) +
  geom_vline(xintercept = 0) +
  geom_vline(aes(xintercept = avg_util), linetype = 2) + 
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), limit = c(0, 200)) + 
  scale_fill_brewer(palette = "Set1") +
  scale_colour_brewer(palette = "Set1") +
  labs(x = "Difference in QALYs", y = "") +
  facet_wrap(vars(strategy), nrow = 1) + 
  theme_linedraw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

psa_plot1 <- psa_plot1.1 + psa_plot1.2 + 
  plot_annotation(tag_levels = 'A')

# Plot 2
psa_plot2.1 <- psa_plot_data |> 
  filter(strategy != "Treatment to Observation Delta") |> 
  ggplot(aes(le)) + 
  geom_histogram(aes(fill = strategy, colour = strategy),
                 bins = 30, alpha = 0.8, position = "identity") +
  geom_vline(aes(xintercept = avg_le), linetype = 2) + 
  scale_fill_brewer(palette = "Set2") +
  scale_colour_brewer(palette = "Set2") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), limit = c(0, 200)) +
  facet_wrap(vars(strategy), nrow = 2) + 
  labs(x = "Life Expectancy (Years)", y = "Count", 
       fill = "Strategy") +
  theme_linedraw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )
  

psa_plot2.2 <- psa_plot_data |>
  filter(strategy == "Treatment to Observation Delta") |> 
  ggplot(aes(le)) + 
  geom_histogram(
    colour = "#9999CC", fill = "#9999CC",
    alpha = 0.8, bins = 30, 
    position = "identity") +
  geom_vline(xintercept = 0) +
  geom_vline(aes(xintercept = avg_le), linetype = 2) + 
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), limit = c(0, 200)) + 
  scale_fill_brewer(palette = "Set1") +
  scale_colour_brewer(palette = "Set1") +
  labs(x = "Difference in Life Expectancy (Years)", y = "") +
  facet_wrap(vars(strategy), nrow = 1) + 
  theme_linedraw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

psa_plot2 <- psa_plot2.1 + psa_plot2.2 + 
  plot_annotation(tag_levels = 'A')

psa_plot_comb <- (psa_plot1.1 + psa_plot1.2) / 
  (psa_plot2.1 + psa_plot2.2) + 
  plot_annotation(tag_levels = "A")

ggsave("psa_plot1.pdf", plot = psa_plot1, path = "Plots/",
       height = 5, width = 10, units = "in")

ggsave("psa_plot2.pdf", plot = psa_plot2, path = "Plots/",
       height = 5, width = 10, units = "in")

ggsave("psa_plot_comb.pdf", plot = psa_plot_comb, path = "Plots/",
       height = 8, width = 7, units = "in")
