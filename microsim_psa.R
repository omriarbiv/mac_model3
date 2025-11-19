library(tidyverse)
library(parallel)
library(patchwork)
source("microsim.R")
source("Functions/gen_psa_params.R")
source("Functions/functions.R")

colour1 <- c("#2D93AD", "#985F99", "#E9AFA3")
colour2 <- c("#A3333D", "#F2D0A4", "#629677")

###########################################################
### Setup

# save_dir <- "/scratch/oarbiv/mac_model3_results/"
save_dir <- "Results/no_discount/"
# files_dir <- "/home/oarbiv/mac_model3/"
files_dir <- ""
n_cores <- detectCores(logical = FALSE)

baseline.vals <- read.csv(paste0(files_dir, "values/baseline.csv"), 
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

life.table <- read.csv(paste0(files_dir, "values/2022_lifetable.csv")) |>
  select("Age", "qx") |>  # qx is the death rate
  rename(age = Age, control_rate = qx) |>
  # Changing age to be just the number
  mutate(age = as.numeric(str_split_i(age, " ", i = 1))) |>
  as.matrix()

cols_names <- c("trt_util", "obs_util", 
                "trt_le", "obs_le",
                paste0("trt_util_y", 1:40),
                paste0("obs_util_y", 1:40))

###########################################################
### Core function

run_psa <- function(iter, save_psa = T,
                    params_base = params_base, 
                    baseline.vals = baseline.vals, 
                    life.table = life.table, ni = ni) {
  
  # Generating random parameters from their distributions
  d_psa_params <- gen_psa_params(n_psa = iter[length(iter)],
                                 pr = 0,
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
  psa_result <- mclapply(1:length(iter), \(x)
                       cbind(psa.t[[x]]$util_total, 
                             psa.o[[x]]$util_total,
                             life_expectancy(psa.t[[x]]$state_matrix),
                             life_expectancy(psa.o[[x]]$state_matrix),
                             yearly_qaly(psa.t[[x]]$util_matrix),
                             yearly_qaly(psa.o[[x]]$util_matrix)),
                       mc.cores = n_cores - 1)
  
  for(i in 1:length(iter)) {
    colnames(psa_result[[i]]) <- cols_names
  }
  
  if(save_psa) {
    saveRDS(
      psa_result, file = paste0(
      save_dir, "psa/psa_n", 
      iter[1], "_", iter[length(iter)], 
      ".RData")
    )
  }
  return(psa_result)
}

###########################################################
### Running PSA

ni <- 1e4
fl <- 40
niter <- 1e4
move_by <- 50
s_seq <- seq(1, niter, by = move_by)

for (i in 1:length(s_seq)) {
  start_time <- Sys.time()
  iter <- s_seq[i]:(s_seq[i] + move_by - 1)
  set.seed(iter[1])
  cat("Running sequence ", iter[1], " to ", iter[length(iter)], "\n")
  p <- run_psa(iter, params_base = params_base, life.table = life.table)
  print(Sys.time() - start_time)
}

###########################################################
### Analyzing PSA

# Find files in Results folder
psa_files <- str_c(save_dir, "psa/", 
                   list.files(paste0(save_dir, "psa/")))

# Chunk PSA files for memory benefit
chunk_size <- 10
chunk_start <- 1
s_seq <- seq(chunk_start, length(psa_files), by = chunk_size)
psa_mean <- tibble()

for(i in 1:length(s_seq)) {
  files_to_download <- psa_files[s_seq[i] : (s_seq[i] + chunk_size - 1)]
  psa_df <- mclapply(files_to_download, readRDS, mc.cores = n_cores - 1) 
  psa_df_tmp <- map(1:length(psa_df), 
                    ~ bind_rows(map(psa_df[[.x]], colMeans))) |> 
    bind_rows() |> 
    mutate(across(ends_with("util"), ~ .x / 12),
           across(ends_with("le"), ~ .x - 70),
           across(matches(".*_util_y\\d.*"), ~ .x / 12))
  
  psa_mean <- bind_rows(psa_mean, psa_df_tmp)
}

psa_mean <- psa_mean |> 
  mutate(
    delta_util = trt_util - obs_util,
    delta_le = trt_le - obs_le,
    delta_util_y1 = trt_util_y1 - obs_util_y1,
    delta_util_y2 = trt_util_y2 - obs_util_y2,
    delta_util_y3 = trt_util_y3 - obs_util_y3,
    delta_util_y4 = trt_util_y4 - obs_util_y4,
    delta_util_y5 = trt_util_y5 - obs_util_y5,
    delta_util_y6 = trt_util_y6 - obs_util_y6,
    delta_util_y7 = trt_util_y7 - obs_util_y7,
    delta_util_y8 = trt_util_y8 - obs_util_y8,
    delta_util_y9 = trt_util_y9 - obs_util_y9,
    delta_util_y10 = trt_util_y10 - obs_util_y10,
    delta_util_y11 = trt_util_y11 - obs_util_y11,
    delta_util_y12 = trt_util_y12 - obs_util_y12,
    delta_util_y13 = trt_util_y13 - obs_util_y13,
    delta_util_y14 = trt_util_y14 - obs_util_y14,
    delta_util_y15 = trt_util_y15 - obs_util_y15,
    delta_util_y16 = trt_util_y16 - obs_util_y16,
    delta_util_y17 = trt_util_y17 - obs_util_y17,
    delta_util_y18 = trt_util_y18 - obs_util_y18,
    delta_util_y19 = trt_util_y19 - obs_util_y19,
    delta_util_y20 = trt_util_y20 - obs_util_y20,
    delta_util_y21 = trt_util_y21 - obs_util_y21,
    delta_util_y22 = trt_util_y22 - obs_util_y22,
    delta_util_y23 = trt_util_y23 - obs_util_y23,
    delta_util_y24 = trt_util_y24 - obs_util_y24,
    delta_util_y25 = trt_util_y25 - obs_util_y25,
    delta_util_y26 = trt_util_y26 - obs_util_y26,
    delta_util_y27 = trt_util_y27 - obs_util_y27,
    delta_util_y28 = trt_util_y28 - obs_util_y28,
    delta_util_y29 = trt_util_y29 - obs_util_y29,
    delta_util_y30 = trt_util_y30 - obs_util_y30,
    delta_util_y31 = trt_util_y31 - obs_util_y31,
    delta_util_y32 = trt_util_y32 - obs_util_y32,
    delta_util_y33 = trt_util_y33 - obs_util_y33,
    delta_util_y34 = trt_util_y34 - obs_util_y34,
    delta_util_y35 = trt_util_y35 - obs_util_y35,
    delta_util_y36 = trt_util_y36 - obs_util_y36,
    delta_util_y37 = trt_util_y37 - obs_util_y37,
    delta_util_y38 = trt_util_y38 - obs_util_y38,
    delta_util_y39 = trt_util_y39 - obs_util_y39,
    delta_util_y40 = trt_util_y40 - obs_util_y40
  )

psa_mean_total <- psa_mean |> 
  select(!matches("y\\d.*$"))

psa_year_mean <- psa_mean |> 
  select(matches("y\\d.*$")) |> 
  summarise(across(everything(), mean)) |> 
  mutate(measure = "mean")

psa_year_lci <- psa_mean |> 
  select(matches("y\\d.*$")) |> 
  summarise(across(everything(), ~ quantile(.x, 0.025))) |> 
  mutate(measure = "lci")

psa_year_uci <- psa_mean |> 
  select(matches("y\\d.*$")) |> 
  summarise(across(everything(), ~ quantile(.x, 0.975))) |> 
  mutate(measure = "uci")

psa_year <- bind_rows(
  psa_year_mean,
  psa_year_lci,
  psa_year_uci
) |> 
  pivot_longer(
    cols = !measure,
    names_to = c("strategy", "year"),
    names_pattern = "(.*)_util_y(.*)"
  ) |> 
  pivot_wider(
    names_from = measure,
    values_from = value
  ) |> 
  add_row(
    strategy = "trt", year = "0", 
    mean = 0, lci = 0, uci = 0,
  ) |> 
  add_row(
    strategy = "obs", year = "0", 
    mean = 0, lci = 0, uci = 0,
  ) |> 
  add_row(
    strategy = "delta", year = "0", 
    mean = 0, lci = 0, uci = 0,
  ) |> 
  mutate(
    year = as.numeric(year),
    Strategy = case_when(
      strategy == "obs" ~ "Observation",
      strategy == "trt" ~ "Treatment",
      strategy == "delta" ~ "Difference"
    )
  ) 


# Read in raw result and obtain means and deltas
# psa_result <- map(psa_files, readRDS) |> 
#   tibble() |> 
#   rename("raw" = 1) |> 
#   unnest(raw) |> 
#   mutate(
#     trt_util = map(raw, ~ mean(.x[, "trt_util"])),
#     obs_util = map(raw, ~ mean(.x[, "obs_util"])),
#     delta_util = map(raw, ~ mean(.x[, "trt_util"] - .x[, "obs_util"])),
#     trt_le = map(raw, ~ mean(.x[, "trt_le"])),
#     obs_le = map(raw, ~ mean(.x[, "obs_le"])),
#     delta_le = map(raw, ~ mean(.x[, "trt_le"] - .x[, "obs_le"])),
#   ) |> 
#   unnest(!1) |> 
#   mutate(across(!raw, ~ .x / 12))

###########################################################
### Analyzing PSA

psa_summary <- psa_mean_total |> 
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
  select(!c(util, le)) |>
  # relocate(starts_with("mcse"), .after = "delta_le") |> 
  slice(1:3) |>
  pivot_longer(
    cols = !1,
    names_to = c(".value", "outcome"),
    names_sep = "_"
  ) |> 
  arrange(desc(outcome), desc(strategy))

# Plot data
psa_plot_data <- psa_mean_total |> 
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
                         strategy == "delta" ~ 
                           "Treatment to Observation Difference"))

psa_plot1.1 <- psa_plot_data |> 
  filter(strategy != "Treatment to Observation Difference") |> 
  ggplot(aes(util)) + 
  geom_histogram(aes(fill = strategy, colour = strategy),
                 bins = 30, alpha = 0.8, position = "identity") +
  geom_vline(aes(xintercept = avg_util), linetype = 2) + 
  scale_colour_manual(values = colour1, aesthetics = c("colour", "fill")) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 25)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 2500)) +
  facet_wrap(vars(strategy), nrow = 2, scales = "free") + 
  labs(x = "QALYs", y = "Count", fill = "Strategy") +
  # theme_linedraw() +
  theme_classic() +
  theme(
    strip.background = element_rect(colour = "white", fill = "white"),
    strip.text = element_text(size = 14),
    axis.line = element_line(),
    # panel.grid.major = element_blank(),
    # panel.grid.minor = element_blank(),
    legend.position = "none",
    axis.ticks.length = unit(5, "points")
  )

psa_plot1.2 <- psa_plot_data |>
  filter(strategy == "Treatment to Observation Difference") |> 
  mutate(strategy = "Treatment to Observation\nDifference") |> 
  # mutate(util0 = ifelse(util > 0, "G", "L")) |> 
  ggplot(aes(util)) + 
  geom_histogram(
    # aes(colour = util0, fill = util0),
    colour = colour1[3], fill = colour1[3],
    alpha = 0.8, bins = 30, 
    position = "identity") +
  # geom_vline(aes(xintercept = avg)) +
  geom_vline(xintercept = 0) +
  geom_vline(aes(xintercept = avg_util), linetype = 2) + 
  scale_x_continuous(expand = c(0, 0), limits = c(-2, 8)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 2500)) +
  scale_fill_brewer(palette = "Set1") +
  scale_colour_brewer(palette = "Set1") +
  labs(x = "Difference in QALYs", y = "") +
  facet_wrap(vars(strategy), nrow = 1) + 
  theme_classic() +
  theme(
    strip.background = element_rect(colour = "white", fill = "white"),
    strip.text = element_text(size = 14),
    legend.position = "none",
    axis.ticks.length = unit(5, "points")
  )

psa_plot1 <- psa_plot1.1 + psa_plot1.2 + 
  plot_annotation(tag_levels = 'A')

# Plot 2
psa_plot2.1 <- psa_plot_data |> 
  filter(strategy != "Treatment to Observation Difference") |> 
  ggplot(aes(le)) + 
  geom_histogram(aes(fill = strategy, colour = strategy),
                 bins = 30, alpha = 0.8, position = "identity") +
  geom_vline(aes(xintercept = avg_le), linetype = 2) + 
  scale_colour_manual(values = colour2, aesthetics = c("colour", "fill")) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 25)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 3000)) +
  facet_wrap(vars(strategy), nrow = 2, scales = "free") + 
  labs(x = "Life Expectancy (Years)", y = "Count", 
       fill = "Strategy") +
  theme_classic() +
  theme(
    strip.background = element_rect(colour = "white", fill = "white"),
    strip.text = element_text(size = 14),
    legend.position = "none",
    axis.ticks.length = unit(5, "points")
  )
  

psa_plot2.2 <- psa_plot_data |>
  filter(strategy == "Treatment to Observation Difference") |> 
  mutate(strategy = "Treatment to Observation\nDifference") |> 
  ggplot(aes(le)) + 
  geom_histogram(
    colour = colour2[3], fill = colour2[3],
    alpha = 0.8, bins = 30, 
    position = "identity") +
  geom_vline(xintercept = 0) +
  geom_vline(aes(xintercept = avg_le), linetype = 2) + 
  scale_x_continuous(expand = c(0, 0), limits = c(-2, 8)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 3000)) +
  scale_fill_brewer(palette = "Set1") +
  scale_colour_brewer(palette = "Set1") +
  labs(x = "Difference in Life Expectancy (Years)", y = "") +
  facet_wrap(vars(strategy), nrow = 1) + 
  theme_classic() +
  theme(
    axis.ticks.length = unit(5, "points"),
    strip.background = element_rect(colour = "white", fill = "white"),
    strip.text = element_text(size = 14)
  )

psa_plot2 <- psa_plot2.1 + psa_plot2.2 + 
  plot_annotation(tag_levels = 'A')

psa_plot_comb <- (psa_plot1.1 + psa_plot1.2) / 
  (psa_plot2.1 + psa_plot2.2) + 
  plot_annotation(tag_levels = "A")

psa_plot3.1 <- psa_year |> 
  filter(Strategy != "Difference") |> 
  ggplot(aes(year, mean, 
             group = Strategy, 
             colour = Strategy)) + 
  geom_line(lwd = 1.3, key_glyph = draw_key_rect) + 
  geom_ribbon(aes(ymin = lci, ymax = uci, 
                  fill = Strategy, colour = NULL),
              alpha = 0.1) + 
  scale_x_continuous(expand = c(0, 0), limits = c(0, 40)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 20)) +
  scale_colour_manual(values = colour1, aesthetics = c("colour", "fill")) +
  labs(x = "Year", y = "QALY") + 
  theme_classic() +
  theme(axis.ticks.length = unit(5, "points"))

psa_plot3.2 <- psa_year |> 
  filter(Strategy == "Difference") |> 
  ggplot(aes(year, mean)) + 
  geom_line(lwd = 1.3, key_glyph = draw_key_rect, 
            colour = colour1[3]) + 
  geom_ribbon(aes(ymin = lci, ymax = uci, colour = NULL),
              alpha = 0.2, fill = colour1[3]) + 
  geom_hline(yintercept = 0) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 40)) +
  scale_y_continuous(expand = c(0, 0), limits = c(-1, 5)) +
  labs(x = "Year", y = "Difference in QALY") + 
  theme_classic() +
  theme(axis.ticks.length = unit(5, "points"))

psa_plot3 <- psa_plot3.1 / psa_plot3.2 + 
  plot_annotation(tag_levels = 'A')

ggsave("psa_plot1.pdf", plot = psa_plot1, 
       path = paste0(save_dir, "plots/"),
       height = 5, width = 10, units = "in")

ggsave("psa_plot2.pdf", plot = psa_plot2, 
       path = paste0(save_dir, "plots/"),
       height = 5, width = 10, units = "in")

ggsave("psa_plot_comb.pdf", plot = psa_plot_comb, 
       path = paste0(save_dir, "plots/"),
       height = 8, width = 7, units = "in")

ggsave("psa_plot3.pdf", plot = psa_plot3,
       path = paste0(save_dir, "plots"),
       height = 5, width = 7, units = "in")

