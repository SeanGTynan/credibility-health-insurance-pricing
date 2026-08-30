# ============================================================
# MSc Actuarial Science Thesis
# Credibility-Based Pricing of Corporate Health Insurance Schemes
# Data: MEPS HC-252 Panel 27 Longitudinal File
# ============================================================

# ------------------------------------------------------------
# 0. Packages and settings
# ------------------------------------------------------------

library(haven)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ranger)
library(xgboost)

method_colours <- c(
  "Credibility" = "#619CFF",
  "Manual Rating" = "#F8766D",
  "Full Experience" = "#00BA38"
)
# Set seed so simulations and train/test splitting can be reproduced
set.seed(183)

# Main assumptions used throughout the analysis
scheme_sizes <- c(10, 50, 100, 500, 1500)
age_profiles <- c("Young", "Middle", "Older")
risk_profiles <- c("Low", "Medium", "High")
n_sims <- 500
k_values <- seq(5, 75, by = 5)
k_values_split <- c(
  seq(5, 100, by = 5),
  seq(125, 500, by = 25),
  seq(600, 2000, by = 100),
  seq(2500, 10000, by = 500)
)
target_loading <- 0.10

# Output folders
if (!dir.exists("outputs")) dir.create("outputs")
if (!dir.exists("outputs/tables")) dir.create("outputs/tables", recursive = TRUE)
if (!dir.exists("outputs/figures")) dir.create("outputs/figures", recursive = TRUE)

# Helper to avoid divide-by-zero prices
floor_price <- function(x) pmax(x, 1)


# ------------------------------------------------------------
# 1. Load data
# ------------------------------------------------------------

meps <- read_dta("data/h252.dta")

raw_data_check <- tibble(
  rows = nrow(meps),
  columns = ncol(meps)
)


# ------------------------------------------------------------
# 2. Clean data
# ------------------------------------------------------------

df <- meps |>
  transmute(
    id = DUPERSID,
    age = as.numeric(AGE1X),
    sex = as.numeric(SEX),
    physical_health = as.numeric(RTHLTH1),
    mental_health = as.numeric(MNHLTH1),
    diabetes = as.numeric(DIABDXY1_M18),
    asthma = as.numeric(ASTHDXY1),
    insurance = as.numeric(INSCOVY1),
    claims_2022 = as.numeric(TOTEXPY1),
    claims_2023 = as.numeric(TOTEXPY2),
    IP_2022 = as.numeric(IPTEXPY1),
    IP_2023 = as.numeric(IPTEXPY2),
    OP_2022 = as.numeric(OPTEXPY1),
    OP_2023 = as.numeric(OPTEXPY2)
  ) |>
  mutate(
    age = ifelse(age < 0, NA, age),
    sex = ifelse(sex < 0, NA, sex),
    physical_health = ifelse(physical_health < 0, NA, physical_health),
    mental_health = ifelse(mental_health < 0, NA, mental_health),
    diabetes = ifelse(diabetes < 0, NA, diabetes),
    asthma = ifelse(asthma < 0, NA, asthma),
    insurance = ifelse(insurance < 0, NA, insurance)
  ) |>
  filter(
    !is.na(age),
    age >= 18,
    !is.na(claims_2022),
    !is.na(claims_2023),
    claims_2022 >= 0,
    claims_2023 >= 0
  )

clean_data_check <- tibble(
  rows = nrow(df),
  columns = ncol(df),
  min_age = min(df$age),
  max_age = max(df$age),
  mean_claims_2022 = mean(df$claims_2022),
  mean_claims_2023 = mean(df$claims_2023)
)


# ------------------------------------------------------------
# 3. Exploratory analysis
# ------------------------------------------------------------

claims_summary <- tibble(
  year = c("2022", "2023"),
  mean_claims = c(mean(df$claims_2022), mean(df$claims_2023)),
  median_claims = c(median(df$claims_2022), median(df$claims_2023)),
  p95_claims = c(quantile(df$claims_2022, 0.95), quantile(df$claims_2023, 0.95)),
  max_claims = c(max(df$claims_2022), max(df$claims_2023))
)

fig_claims_distribution <- ggplot(df, aes(x = claims_2022 + 1)) +
  geom_histogram(
    bins = 50,
    fill = "#619CFF",
    colour = "white",
    linewidth = 0.15
  ) +
  scale_x_log10(
    breaks = c(1, 10, 100, 1000, 10000, 100000, 1000000),
    labels = scales::label_comma()
  ) +
  labs(
    subtitle = "Claims are shown on a logarithmic scale.",
    x = "2022 individual claims",
    y = "Number of individuals"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.subtitle = element_text(size = 10),
    panel.grid.minor = element_blank()
  )


# ------------------------------------------------------------
# 4. Risk score model
# ------------------------------------------------------------

df_model <- df |>
  mutate(
    claims_2022_log = log(claims_2022 + 1),
    sex = factor(sex),
    physical_health = factor(physical_health),
    mental_health = factor(mental_health),
    diabetes = factor(diabetes),
    asthma = factor(asthma),
    insurance = factor(insurance)
  ) |>
  filter(
    !is.na(age),
    !is.na(sex),
    !is.na(physical_health),
    !is.na(mental_health),
    !is.na(diabetes),
    !is.na(asthma),
    !is.na(insurance),
    !is.na(claims_2022_log)
  )

risk_model <- lm(
  claims_2022_log ~ age + sex + physical_health + mental_health + diabetes + asthma + insurance,
  data = df_model
)

df_model <- df_model |>
  mutate(
    risk_score_model = exp(predict(risk_model)) - 1
  )

df_clean <- df_model |>
  mutate(
    age_profile = case_when(
      age < 35 ~ "Young",
      age >= 35 & age < 55 ~ "Middle",
      age >= 55 ~ "Older"
    )
  ) |>
  group_by(age_profile) |>
  mutate(
    risk_band_model = ntile(risk_score_model, 3),
    risk_profile_model = case_when(
      risk_band_model == 1 ~ "Low",
      risk_band_model == 2 ~ "Medium",
      TRUE ~ "High"
    )
  ) |>
  ungroup() |>
  mutate(
    age_profile = factor(age_profile, levels = c("Young", "Middle", "Older")),
    risk_profile_model = factor(risk_profile_model, levels = c("Low", "Medium", "High"))
  )

risk_band_check <- df_clean |>
  group_by(age_profile, risk_profile_model) |>
  summarise(
    members = n(),
    avg_risk_score = mean(risk_score_model),
    avg_claims_2022 = mean(claims_2022),
    avg_claims_2023 = mean(claims_2023),
    .groups = "drop"
  )

# Risk profile cutoffs by age group
risk_cutoffs <- df_clean |>
  group_by(age_profile) |>
  summarise(
    low_medium_cutoff = quantile(risk_score_model, 1/3, na.rm = TRUE),
    medium_high_cutoff = quantile(risk_score_model, 2/3, na.rm = TRUE),
    min_score = min(risk_score_model, na.rm = TRUE),
    max_score = max(risk_score_model, na.rm = TRUE),
    .groups = "drop"
  )


# ------------------------------------------------------------
# 5. Scheme simulation
# ------------------------------------------------------------

simulate_scheme <- function(size, age_prof, risk_prof, sim_id) {
  pool <- df_clean |>
    filter(
      age_profile == age_prof,
      risk_profile_model == risk_prof
    )
  
  sample_members <- pool |>
    slice_sample(n = size, replace = TRUE)
  
  tibble(
    sim_id = sim_id,
    scheme_size = size,
    age_profile = age_prof,
    risk_profile = risk_prof,
    claims_2022_pm = mean(sample_members$claims_2022),
    claims_2023_pm = mean(sample_members$claims_2023),
    IP_2022_pm = mean(sample_members$IP_2022),
    IP_2023_pm = mean(sample_members$IP_2023),
    OP_2022_pm = mean(sample_members$OP_2022),
    OP_2023_pm = mean(sample_members$OP_2023),
    total_claims_2022 = sum(sample_members$claims_2022),
    total_claims_2023 = sum(sample_members$claims_2023)
  )
}

sim_results <- expand.grid(
  scheme_size = scheme_sizes,
  age_profile = age_profiles,
  risk_profile = risk_profiles,
  sim_id = 1:n_sims
) |>
  rowwise() |>
  do(simulate_scheme(
    size = .$scheme_size,
    age_prof = .$age_profile,
    risk_prof = .$risk_profile,
    sim_id = .$sim_id
  )) |>
  ungroup()

scheme_summary <- sim_results |>
  group_by(scheme_size, age_profile, risk_profile) |>
  summarise(
    avg_2022_claims = mean(claims_2022_pm, na.rm = TRUE),
    avg_2023_claims = mean(claims_2023_pm, na.rm = TRUE),
    sd_2023_claims = sd(claims_2023_pm, na.rm = TRUE),
    p10_2023_claims = quantile(claims_2023_pm, 0.10, na.rm = TRUE),
    median_2023_claims = median(claims_2023_pm, na.rm = TRUE),
    p90_2023_claims = quantile(claims_2023_pm, 0.90, na.rm = TRUE),
    avg_IP_2023 = mean(IP_2023_pm, na.rm = TRUE),
    avg_OP_2023 = mean(OP_2023_pm, na.rm = TRUE),
    .groups = "drop"
  )

scheme_size_volatility <- sim_results |>
  group_by(scheme_size) |>
  summarise(
    sd_2023_claims = sd(claims_2023_pm, na.rm = TRUE),
    avg_2023_claims = mean(claims_2023_pm, na.rm = TRUE),
    .groups = "drop"
  )

fig_scheme_volatility <- ggplot(
  scheme_size_volatility,
  aes(x = factor(scheme_size), y = sd_2023_claims)
) +
  geom_col(fill = "#619CFF", width = 0.65) +
  geom_text(aes(label = round(sd_2023_claims, 0)),
            vjust = -0.4, size = 4) +
  labs(
    subtitle = "Volatility is measured using the standard deviation\nof simulated 2023 claims per member.",
    x = "Scheme size",
    y = "Volatility of 2023 claims per member"
  ) +
  scale_y_continuous(labels = scales::comma) +
  theme_minimal(base_size = 14) +
  theme(
    plot.subtitle = element_text(size = 11),
    axis.title = element_text(face = "bold"),
    plot.margin = margin(10, 20, 10, 10)
  ) +
  expand_limits(y = max(scheme_size_volatility$sd_2023_claims) * 1.08)


scheme_claim_ranges <- sim_results |>
  group_by(scheme_size) |>
  summarise(
    p10_2023_claims = quantile(claims_2023_pm, 0.10, na.rm = TRUE),
    median_2023_claims = median(claims_2023_pm, na.rm = TRUE),
    p90_2023_claims = quantile(claims_2023_pm, 0.90, na.rm = TRUE),
    .groups = "drop"
  )

fig_scheme_claim_range <- ggplot(
  scheme_claim_ranges,
  aes(x = factor(scheme_size))
) +
  geom_linerange(
    aes(ymin = p10_2023_claims, ymax = p90_2023_claims),
    colour = "#9EC5FE",
    linewidth = 1.4
  ) +
  geom_point(
    aes(y = median_2023_claims),
    colour = "#619CFF",
    size = 3.5
  ) +
  geom_text(
    aes(
      y = median_2023_claims,
      label = scales::comma(round(median_2023_claims, 0))
    ),
    vjust = -1,
    size = 3
  ) +
  scale_y_continuous(labels = scales::comma) +
  labs(
    subtitle = "Lines show the 10th to 90th percentile. Dots show the median outcome.",
    x = "Scheme size",
    y = "2023 claims per member"
  ) +
  theme_minimal(base_size = 12)

# ------------------------------------------------------------
# 6. Pricing methods
# ------------------------------------------------------------

market_2022 <- mean(df_clean$claims_2022)
market_2023 <- mean(df_clean$claims_2023)
trend_factor <- market_2023 / market_2022

manual_rates <- df_clean |>
  group_by(age_profile, risk_profile_model) |>
  summarise(
    manual_2022 = mean(claims_2022),
    .groups = "drop"
  ) |>
  rename(risk_profile = risk_profile_model)

apply_total_pricing <- function(k_value) {
  sim_results |>
    left_join(manual_rates, by = c("age_profile", "risk_profile")) |>
    mutate(
      k = k_value,
      Z = scheme_size / (scheme_size + k),
      price_manual = floor_price(manual_2022 * trend_factor),
      price_experience = floor_price(claims_2022_pm * trend_factor),
      price_credibility = floor_price((Z * claims_2022_pm + (1 - Z) * manual_2022) * trend_factor),
      error_manual = price_manual - claims_2023_pm,
      error_experience = price_experience - claims_2023_pm,
      error_credibility = price_credibility - claims_2023_pm,
      abs_error_manual = abs(error_manual),
      abs_error_experience = abs(error_experience),
      abs_error_credibility = abs(error_credibility),
      lr_manual = claims_2023_pm / price_manual,
      lr_experience = claims_2023_pm / price_experience,
      lr_credibility = claims_2023_pm / price_credibility
    )
}

make_pricing_long <- function(pricing_data) {
  pricing_data |>
    select(
      scheme_size, age_profile, risk_profile, sim_id,
      abs_error_manual, abs_error_experience, abs_error_credibility
    ) |>
    pivot_longer(
      cols = c(abs_error_manual, abs_error_experience, abs_error_credibility),
      names_to = "method",
      values_to = "abs_error"
    ) |>
    mutate(
      method = case_when(
        method == "abs_error_manual" ~ "Manual Rating",
        method == "abs_error_experience" ~ "Full Experience",
        method == "abs_error_credibility" ~ "Credibility"
      )
    )
}


# ------------------------------------------------------------
# 7. Credibility sensitivity
# ------------------------------------------------------------

evaluate_credibility_k <- function(k_value) {
  apply_total_pricing(k_value) |>
    group_by(scheme_size) |>
    summarise(
      k = k_value,
      mean_abs_error_credibility = mean(abs_error_credibility),
      median_abs_error_credibility = median(abs_error_credibility),
      .groups = "drop"
    )
}

cred_sensitivity <- bind_rows(lapply(k_values, evaluate_credibility_k))

cred_sensitivity_overall <- cred_sensitivity |>
  group_by(k) |>
  summarise(
    overall_mean_abs_error = mean(mean_abs_error_credibility),
    overall_median_abs_error = mean(median_abs_error_credibility),
    .groups = "drop"
  )

best_total_k <- cred_sensitivity_overall |>
  arrange(overall_mean_abs_error) |>
  slice(1)

k_final_total <- best_total_k$k
pricing_results <- apply_total_pricing(k_final_total)
pricing_long <- make_pricing_long(pricing_results)

method_summary <- pricing_long |>
  group_by(scheme_size, method) |>
  summarise(
    mean_abs_error = mean(abs_error),
    median_abs_error = median(abs_error),
    sd_abs_error = sd(abs_error),
    p5_abs_error = quantile(abs_error, 0.05),
    p95_abs_error = quantile(abs_error, 0.95),
    .groups = "drop"
  )

loss_ratio_summary <- pricing_results |>
  group_by(scheme_size, age_profile, risk_profile)|>
  summarise(
    avg_LR_manual = mean(lr_manual),
    avg_LR_experience = mean(lr_experience),
    avg_LR_credibility = mean(lr_credibility),
    sd_LR_manual = sd(lr_manual),
    sd_LR_experience = sd(lr_experience),
    sd_LR_credibility = sd(lr_credibility),
    .groups = "drop"
  )

fig_cred_sensitivity <- ggplot(cred_sensitivity, aes(x = k, y = mean_abs_error_credibility)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ scheme_size, scales = "free_y") +
  labs(
    title = "Sensitivity of Credibility Pricing Accuracy to Credibility Parameter",
    x = "Credibility parameter k",
    y = "Mean absolute pricing error"
  )

fig_pricing_accuracy <- ggplot(method_summary, aes(x = factor(scheme_size), y = mean_abs_error, fill = method)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = method_colours) +
  labs(
    x = "Scheme size",
    y = "Mean absolute pricing error per member",
    fill = "Pricing method"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

fig_pricing_profile <- pricing_long |>
  group_by(scheme_size, age_profile, risk_profile, method) |>
  summarise(mean_abs_error = mean(abs_error), .groups = "drop") |>
  ggplot(aes(x = risk_profile, y = mean_abs_error, fill = method)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = method_colours) +
  facet_grid(age_profile ~ scheme_size) +
  labs(
    title = "Pricing Accuracy by Age Profile, Risk Profile and Scheme Size",
    x = "Risk profile",
    y = "Mean absolute pricing error",
    fill = "Pricing method"
  )

# ------------------------------------------------------------
# 8. Sensitivity: illustrative mixed-composition adult schemes
# ------------------------------------------------------------

mixed_profiles <- tibble::tribble(
  ~profile_name,        ~age_profile, ~risk_profile, ~share,
  
  "Young workforce",    "Young",      "Low",         0.35,
  "Young workforce",    "Young",      "Medium",      0.30,
  "Young workforce",    "Young",      "High",        0.10,
  "Young workforce",    "Middle",      "Low",         0.10,
  "Young workforce",    "Middle",      "Medium",      0.10,
  "Young workforce",    "Older",      "Medium",      0.05,
  
  "Balanced workforce", "Young",      "Low",         0.15,
  "Balanced workforce", "Young",      "Medium",      0.15,
  "Balanced workforce", "Middle",      "Low",         0.15,
  "Balanced workforce", "Middle",      "Medium",      0.25,
  "Balanced workforce", "Middle",      "High",        0.10,
  "Balanced workforce", "Older",      "Medium",      0.15,
  "Balanced workforce", "Older",      "High",        0.05,
  
  "Older workforce",    "Young",      "Medium",      0.05,
  "Older workforce",    "Middle",      "Medium",      0.20,
  "Older workforce",    "Middle",      "High",        0.15,
  "Older workforce",    "Older",      "Low",         0.10,
  "Older workforce",    "Older",      "Medium",      0.25,
  "Older workforce",    "Older",      "High",        0.25
)

# Check shares sum to 1
mixed_profile_check <- mixed_profiles |>
  group_by(profile_name) |>
  summarise(total_share = sum(share), .groups = "drop") 

simulate_mixed_scheme <- function(size, profile, sim_id) {
  
  profile_weights <- mixed_profiles |>
    filter(profile_name == profile) |>
    mutate(
      raw_n = size * share,
      n_members = floor(raw_n),
      remainder = raw_n - n_members
    )
  
  # Adjust rounding so total sampled members equals scheme size
  shortfall <- size - sum(profile_weights$n_members)
  
  if (shortfall > 0) {
    add_rows <- order(profile_weights$remainder, decreasing = TRUE)[1:shortfall]
    profile_weights$n_members[add_rows] <- profile_weights$n_members[add_rows] + 1
  }
  
  sampled_members <- profile_weights |>
    filter(n_members > 0) |>
    group_by(age_profile, risk_profile, n_members) |>
    group_modify(~ {
      df_clean |>
        filter(
          age_profile == .y$age_profile,
          risk_profile_model == .y$risk_profile
        ) |>
        slice_sample(n = .y$n_members, replace = TRUE) |>
        select(-age_profile)
    }) |>
    ungroup()
  
  tibble(
    sim_id = sim_id,
    scheme_size = size,
    profile_name = profile,
    
    claims_2022_pm = mean(sampled_members$claims_2022),
    claims_2023_pm = mean(sampled_members$claims_2023),
    
    IP_2022_pm = mean(sampled_members$IP_2022),
    IP_2023_pm = mean(sampled_members$IP_2023),
    
    OP_2022_pm = mean(sampled_members$OP_2022),
    OP_2023_pm = mean(sampled_members$OP_2023)
  )
}

set.seed(183)

mixed_sim_results <- expand.grid(
  scheme_size = scheme_sizes,
  profile_name = unique(mixed_profiles$profile_name),
  sim_id = 1:n_sims
) |>
  rowwise() |>
  do(simulate_mixed_scheme(
    size = .$scheme_size,
    profile = .$profile_name,
    sim_id = .$sim_id
  )) |>
  ungroup()

# Manual rate for each mixed profile = weighted average of age/risk manual rates

mixed_manual_rates <- mixed_profiles |>
  left_join(manual_rates, by = c("age_profile", "risk_profile")) |>
  group_by(profile_name) |>
  summarise(
    mixed_manual_2022 = sum(share * manual_2022),
    .groups = "drop"
  )

mixed_pricing_results <- mixed_sim_results |>
  left_join(mixed_manual_rates, by = "profile_name") |>
  mutate(
    Z = scheme_size / (scheme_size + k_final_total),
    
    price_manual = mixed_manual_2022 * trend_factor,
    price_experience = claims_2022_pm * trend_factor,
    price_credibility = (Z * claims_2022_pm + (1 - Z) * mixed_manual_2022) * trend_factor,
    
    price_manual = floor_price(price_manual),
    price_experience = floor_price(price_experience),
    price_credibility = floor_price(price_credibility),
    
    abs_error_manual = abs(price_manual - claims_2023_pm),
    abs_error_experience = abs(price_experience - claims_2023_pm),
    abs_error_credibility = abs(price_credibility - claims_2023_pm),
    
    lr_manual = claims_2023_pm / price_manual,
    lr_experience = claims_2023_pm / price_experience,
    lr_credibility = claims_2023_pm / price_credibility
  )

mixed_method_summary <- mixed_pricing_results |>
  group_by(scheme_size, profile_name) |>
  summarise(
    Manual_Rating = mean(abs_error_manual),
    Full_Experience = mean(abs_error_experience),
    Credibility = mean(abs_error_credibility),
    .groups = "drop"
  )

mixed_long <- mixed_method_summary |>
  pivot_longer(
    cols = c(Manual_Rating, Full_Experience, Credibility),
    names_to = "method",
    values_to = "mean_abs_error"
  ) |>
  mutate(
    method = case_when(
      method == "Manual_Rating" ~ "Manual Rating",
      method == "Full_Experience" ~ "Full Experience",
      method == "Credibility" ~ "Credibility"
    ),
    method = factor(
      method,
      levels = c("Manual Rating", "Full Experience", "Credibility")
    ),
    profile_name = factor(
      profile_name,
      levels = c("Young workforce", "Balanced workforce", "Older workforce")
    )
  )

fig_mixed_pricing_accuracy <- ggplot(
  mixed_long,
  aes(x = factor(scheme_size), y = mean_abs_error, fill = method)
) +
  geom_col(position = "dodge", width = 0.75) +
  scale_fill_manual(values = method_colours) +
  facet_wrap(~ profile_name, nrow = 1) +
  scale_y_continuous(labels = scales::comma) +
  labs(
    subtitle = "Lower values indicate better pricing accuracy.",
    x = "Scheme size",
    y = "Mean absolute pricing error per member",
    fill = "Pricing method"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold")
  )


mixed_method_winners <- mixed_method_summary |>
  rowwise() |>
  mutate(
    best_method = names(which.min(c(
      Manual_Rating = Manual_Rating,
      Full_Experience = Full_Experience,
      Credibility = Credibility
    ))),
    best_error = min(c(Manual_Rating, Full_Experience, Credibility)),
    credibility_vs_manual_improvement =
      (Manual_Rating - Credibility) / Manual_Rating
  ) |>
  ungroup()


mixed_winner_count <- mixed_method_winners |>
  count(best_method)


# ------------------------------------------------------------
# 8. IP/OP analysis
# ------------------------------------------------------------

IP_trend <- mean(df_clean$IP_2023) / mean(df_clean$IP_2022)
OP_trend <- mean(df_clean$OP_2023) / mean(df_clean$OP_2022)

manual_rates_split <- df_clean |>
  group_by(age_profile, risk_profile_model) |>
  summarise(
    manual_IP_2022 = mean(IP_2022),
    manual_OP_2022 = mean(OP_2022),
    .groups = "drop"
  ) |>
  rename(risk_profile = risk_profile_model)

pricing_split <- sim_results |>
  left_join(manual_rates_split, by = c("age_profile", "risk_profile")) |>
  mutate(
    Z = scheme_size / (scheme_size + k_final_total),
    price_IP_manual = manual_IP_2022 * IP_trend,
    price_IP_experience = IP_2022_pm * IP_trend,
    price_IP_credibility = (Z * IP_2022_pm + (1 - Z) * manual_IP_2022) * IP_trend,
    price_OP_manual = manual_OP_2022 * OP_trend,
    price_OP_experience = OP_2022_pm * OP_trend,
    price_OP_credibility = (Z * OP_2022_pm + (1 - Z) * manual_OP_2022) * OP_trend,
    abs_error_IP_manual = abs(price_IP_manual - IP_2023_pm),
    abs_error_IP_experience = abs(price_IP_experience - IP_2023_pm),
    abs_error_IP_credibility = abs(price_IP_credibility - IP_2023_pm),
    abs_error_OP_manual = abs(price_OP_manual - OP_2023_pm),
    abs_error_OP_experience = abs(price_OP_experience - OP_2023_pm),
    abs_error_OP_credibility = abs(price_OP_credibility - OP_2023_pm)
  )

split_summary <- pricing_split |>
  group_by(scheme_size) |>
  summarise(
    IP_manual_MAE = mean(abs_error_IP_manual),
    IP_experience_MAE = mean(abs_error_IP_experience),
    IP_credibility_MAE = mean(abs_error_IP_credibility),
    OP_manual_MAE = mean(abs_error_OP_manual),
    OP_experience_MAE = mean(abs_error_OP_experience),
    OP_credibility_MAE = mean(abs_error_OP_credibility),
    .groups = "drop"
  )

split_long <- pricing_split |>
  select(
    scheme_size, age_profile, risk_profile, sim_id,
    abs_error_IP_manual, abs_error_IP_experience, abs_error_IP_credibility,
    abs_error_OP_manual, abs_error_OP_experience, abs_error_OP_credibility
  ) |>
  pivot_longer(
    cols = starts_with("abs_error"),
    names_to = "metric",
    values_to = "abs_error"
  ) |>
  mutate(
    claim_type = case_when(
      grepl("_IP_", metric) ~ "Inpatient",
      grepl("_OP_", metric) ~ "Outpatient"
    ),
    method = case_when(
      grepl("manual", metric) ~ "Manual Rating",
      grepl("experience", metric) ~ "Full Experience",
      grepl("credibility", metric) ~ "Credibility"
    )
  )

split_method_summary <- split_long |>
  group_by(scheme_size, claim_type, method) |>
  summarise(
    mean_abs_error = mean(abs_error),
    median_abs_error = median(abs_error),
    sd_abs_error = sd(abs_error),
    .groups = "drop"
  )

evaluate_split_k <- function(k_value) {
  sim_results |>
    left_join(manual_rates_split, by = c("age_profile", "risk_profile")) |>
    mutate(
      k = k_value,
      Z = scheme_size / (scheme_size + k),
      price_IP_credibility = (Z * IP_2022_pm + (1 - Z) * manual_IP_2022) * IP_trend,
      price_OP_credibility = (Z * OP_2022_pm + (1 - Z) * manual_OP_2022) * OP_trend,
      abs_error_IP_credibility = abs(price_IP_credibility - IP_2023_pm),
      abs_error_OP_credibility = abs(price_OP_credibility - OP_2023_pm)
    ) |>
    group_by(scheme_size) |>
    summarise(
      k = k_value,
      IP_MAE = mean(abs_error_IP_credibility),
      OP_MAE = mean(abs_error_OP_credibility),
      .groups = "drop"
    )
}

split_k_sensitivity <- bind_rows(lapply(k_values_split, evaluate_split_k))

best_split_k <- split_k_sensitivity |>
  group_by(scheme_size) |>
  summarise(
    best_IP_k = k[which.min(IP_MAE)],
    best_OP_k = k[which.min(OP_MAE)],
    best_IP_MAE = min(IP_MAE),
    best_OP_MAE = min(OP_MAE),
    .groups = "drop"
  )

split_k_long <- split_k_sensitivity |>
  pivot_longer(
    cols = c(IP_MAE, OP_MAE),
    names_to = "claim_type",
    values_to = "mean_abs_error"
  ) |>
  mutate(
    claim_type = case_when(
      claim_type == "IP_MAE" ~ "Inpatient",
      claim_type == "OP_MAE" ~ "Outpatient"
    )
  )

fig_split_accuracy <- ggplot(
  split_method_summary,
  aes(x = factor(scheme_size), y = mean_abs_error, fill = method)
) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = method_colours) +
  facet_wrap(~ claim_type, scales = "free_y") +
  labs(
    x = "Scheme size",
    y = "Mean absolute pricing error per member",
    fill = "Pricing method"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

fig_split_k <- ggplot(split_k_long, aes(x = k, y = mean_abs_error)) +
  geom_line() +
  geom_point(size = 1.3) +
  scale_x_log10(
    breaks = c(10, 100, 1000, 10000),
    labels = scales::comma
  ) +
  facet_grid(claim_type ~ scheme_size, scales = "free_y") +
  labs(
    subtitle = "Lower values indicate better pricing accuracy.",
    x = "Credibility parameter k, log scale",
    y = "Mean absolute pricing error per member"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold")
  )

# ------------------------------------------------------------
# 9. ML benchmark
# ------------------------------------------------------------

set.seed(183)

model_compare_df <- df_clean |>
  mutate(
    claims_2022_log = log(claims_2022 + 1),
    claims_2023_log = log(claims_2023 + 1),
    has_claim_2023 = ifelse(claims_2023 > 0, 1, 0)
  ) |>
  select(
    claims_2023, claims_2023_log, has_claim_2023,
    claims_2022_log, age, sex, physical_health, mental_health,
    diabetes, asthma, insurance
  ) |>
  na.omit()

train_index <- sample(seq_len(nrow(model_compare_df)), size = 0.7 * nrow(model_compare_df))
train_data <- model_compare_df[train_index, ]
test_data <- model_compare_df[-train_index, ]

model_loglinear <- lm(
  claims_2023_log ~ claims_2022_log + age + sex + physical_health +
    mental_health + diabetes + asthma + insurance,
  data = train_data
)

pred_loglinear <- pmax(exp(predict(model_loglinear, newdata = test_data)) - 1, 0)

model_claim_prob <- glm(
  has_claim_2023 ~ claims_2022_log + age + sex + physical_health +
    mental_health + diabetes + asthma + insurance,
  data = train_data,
  family = binomial()
)

positive_train <- train_data |>
  filter(claims_2023 > 0)

model_claim_amount <- lm(
  log(claims_2023) ~ claims_2022_log + age + sex + physical_health +
    mental_health + diabetes + asthma + insurance,
  data = positive_train
)

prob_claim <- predict(model_claim_prob, newdata = test_data, type = "response")
amount_given_claim <- exp(predict(model_claim_amount, newdata = test_data))
pred_twopart <- pmax(prob_claim * amount_given_claim, 0)

model_rf <- ranger(
  claims_2023_log ~ claims_2022_log + age + sex + physical_health +
    mental_health + diabetes + asthma + insurance,
  data = train_data,
  num.trees = 500,
  mtry = 3,
  min.node.size = 20,
  importance = "impurity",
  seed = 183
)

pred_rf <- pmax(exp(predict(model_rf, data = test_data)$predictions) - 1, 0)

x_train <- model.matrix(
  claims_2023_log ~ claims_2022_log + age + sex + physical_health +
    mental_health + diabetes + asthma + insurance - 1,
  data = train_data
)

x_test <- model.matrix(
  claims_2023_log ~ claims_2022_log + age + sex + physical_health +
    mental_health + diabetes + asthma + insurance - 1,
  data = test_data
)

dtrain <- xgb.DMatrix(data = x_train, label = train_data$claims_2023_log)
dtest <- xgb.DMatrix(data = x_test)

xgb_params <- list(
  objective = "reg:squarederror",
  eval_metric = "mae",
  max_depth = 3,
  eta = 0.05,
  subsample = 0.8,
  colsample_bytree = 0.8
)

xgb_cv <- xgb.cv(
  params = xgb_params,
  data = dtrain,
  nrounds = 500,
  nfold = 5,
  early_stopping_rounds = 25,
  verbose = 0
)

best_nrounds <- xgb_cv$best_iteration

if (is.null(best_nrounds) || length(best_nrounds) == 0 || is.na(best_nrounds)) {
  best_nrounds <- which.min(xgb_cv$evaluation_log$test_mae_mean)
}

model_xgb <- xgb.train(
  params = xgb_params,
  data = dtrain,
  nrounds = best_nrounds,
  verbose = 0
)

pred_xgb <- pmax(exp(predict(model_xgb, dtest)) - 1, 0)

model_comparison <- tibble(
  model = c("Log-linear", "Two-part", "Random Forest", "XGBoost"),
  MAE = c(
    mean(abs(pred_loglinear - test_data$claims_2023)),
    mean(abs(pred_twopart - test_data$claims_2023)),
    mean(abs(pred_rf - test_data$claims_2023)),
    mean(abs(pred_xgb - test_data$claims_2023))
  ),
  RMSE = c(
    sqrt(mean((pred_loglinear - test_data$claims_2023)^2)),
    sqrt(mean((pred_twopart - test_data$claims_2023)^2)),
    sqrt(mean((pred_rf - test_data$claims_2023)^2)),
    sqrt(mean((pred_xgb - test_data$claims_2023)^2))
  )
)

rf_importance <- tibble(
  variable = names(model_rf$variable.importance),
  importance = as.numeric(model_rf$variable.importance)
) |>
  arrange(desc(importance))

xgb_importance <- xgb.importance(
  feature_names = colnames(x_train),
  model = model_xgb
)

# Repeated cross-validation tuning for RF and XGBoost
k_folds <- 5
n_repeats <- 3

rf_grid <- expand.grid(
  mtry = c(4, 5, 6, 7, 8),
  min.node.size = c(50, 100, 150, 200)
)

rf_cv_results <- bind_rows(lapply(seq_len(nrow(rf_grid)), function(i) {
  all_fold_errors <- c()
  
  for (rep_id in 1:n_repeats) {
    fold_id <- sample(rep(1:k_folds, length.out = nrow(train_data)))
    
    for (fold in 1:k_folds) {
      cv_train <- train_data[fold_id != fold, ]
      cv_valid <- train_data[fold_id == fold, ]
      
      rf_temp <- ranger(
        claims_2023_log ~ claims_2022_log + age + sex + physical_health +
          mental_health + diabetes + asthma + insurance,
        data = cv_train,
        num.trees = 300,
        mtry = rf_grid$mtry[i],
        min.node.size = rf_grid$min.node.size[i],
        seed = 183
      )
      
      pred_valid <- pmax(exp(predict(rf_temp, data = cv_valid)$predictions) - 1, 0)
      all_fold_errors <- c(all_fold_errors, mean(abs(pred_valid - cv_valid$claims_2023)))
    }
  }
  
  tibble(
    mtry = rf_grid$mtry[i],
    min.node.size = rf_grid$min.node.size[i],
    repeated_cv_MAE = mean(all_fold_errors),
    repeated_cv_MAE_sd = sd(all_fold_errors)
  )
}))

best_rf <- rf_cv_results |>
  arrange(repeated_cv_MAE) |>
  slice(1)

model_rf_tuned <- ranger(
  claims_2023_log ~ claims_2022_log + age + sex + physical_health +
    mental_health + diabetes + asthma + insurance,
  data = train_data,
  num.trees = 500,
  mtry = best_rf$mtry,
  min.node.size = best_rf$min.node.size,
  importance = "impurity",
  seed = 183
)

pred_rf_tuned <- pmax(exp(predict(model_rf_tuned, data = test_data)$predictions) - 1, 0)

xgb_grid <- expand.grid(
  max_depth = c(1, 2, 3, 4, 5),
  eta = c(0.01, 0.03, 0.05, 0.1),
  subsample = c(0.8),
  colsample_bytree = c(0.8)
)

x_train_full <- model.matrix(
  claims_2023_log ~ claims_2022_log + age + sex + physical_health +
    mental_health + diabetes + asthma + insurance - 1,
  data = train_data
)

y_train_full <- train_data$claims_2023_log
actual_train_full <- train_data$claims_2023

xgb_cv_results <- bind_rows(lapply(seq_len(nrow(xgb_grid)), function(i) {
  all_fold_errors <- c()
  best_rounds <- c()
  
  for (rep_id in 1:n_repeats) {
    fold_id <- sample(rep(1:k_folds, length.out = nrow(train_data)))
    
    for (fold in 1:k_folds) {
      train_rows <- which(fold_id != fold)
      valid_rows <- which(fold_id == fold)
      
      dtrain_cv <- xgb.DMatrix(data = x_train_full[train_rows, ], label = y_train_full[train_rows])
      dvalid_cv <- xgb.DMatrix(data = x_train_full[valid_rows, ], label = y_train_full[valid_rows])
      
      params_temp <- list(
        objective = "reg:squarederror",
        eval_metric = "mae",
        max_depth = xgb_grid$max_depth[i],
        eta = xgb_grid$eta[i],
        subsample = xgb_grid$subsample[i],
        colsample_bytree = xgb_grid$colsample_bytree[i]
      )
      
      xgb_temp <- xgb.train(
        params = params_temp,
        data = dtrain_cv,
        nrounds = 1000,
        watchlist = list(train = dtrain_cv, valid = dvalid_cv),
        early_stopping_rounds = 40,
        verbose = 0
      )
      
      pred_valid <- pmax(exp(predict(xgb_temp, dvalid_cv)) - 1, 0)
      all_fold_errors <- c(all_fold_errors, mean(abs(pred_valid - actual_train_full[valid_rows])))
      
      best_round <- xgb_temp$best_iteration
      if (is.null(best_round) || length(best_round) == 0 || is.na(best_round)) {
        best_round <- 1000
      }
      best_rounds <- c(best_rounds, best_round)
    }
  }
  
  tibble(
    max_depth = xgb_grid$max_depth[i],
    eta = xgb_grid$eta[i],
    subsample = xgb_grid$subsample[i],
    colsample_bytree = xgb_grid$colsample_bytree[i],
    repeated_cv_MAE = mean(all_fold_errors),
    repeated_cv_MAE_sd = sd(all_fold_errors),
    avg_best_round = round(mean(best_rounds, na.rm = TRUE))
  )
}))

best_xgb <- xgb_cv_results |>
  arrange(repeated_cv_MAE) |>
  slice(1)

xgb_params_tuned <- list(
  objective = "reg:squarederror",
  eval_metric = "mae",
  max_depth = best_xgb$max_depth,
  eta = best_xgb$eta,
  subsample = best_xgb$subsample,
  colsample_bytree = best_xgb$colsample_bytree
)

dtrain_full <- xgb.DMatrix(data = x_train_full, label = y_train_full)

model_xgb_tuned <- xgb.train(
  params = xgb_params_tuned,
  data = dtrain_full,
  nrounds = best_xgb$avg_best_round,
  verbose = 0
)

pred_xgb_tuned <- pmax(exp(predict(model_xgb_tuned, dtest)) - 1, 0)

model_comparison_tuned <- tibble(
  model = c("Log-linear", "Two-part", "Random Forest", "Random Forest Tuned", "XGBoost", "XGBoost Tuned"),
  MAE = c(
    mean(abs(pred_loglinear - test_data$claims_2023)),
    mean(abs(pred_twopart - test_data$claims_2023)),
    mean(abs(pred_rf - test_data$claims_2023)),
    mean(abs(pred_rf_tuned - test_data$claims_2023)),
    mean(abs(pred_xgb - test_data$claims_2023)),
    mean(abs(pred_xgb_tuned - test_data$claims_2023))
  ),
  RMSE = c(
    sqrt(mean((pred_loglinear - test_data$claims_2023)^2)),
    sqrt(mean((pred_twopart - test_data$claims_2023)^2)),
    sqrt(mean((pred_rf - test_data$claims_2023)^2)),
    sqrt(mean((pred_rf_tuned - test_data$claims_2023)^2)),
    sqrt(mean((pred_xgb - test_data$claims_2023)^2)),
    sqrt(mean((pred_xgb_tuned - test_data$claims_2023)^2))
  )
)

best_risk_model <- model_comparison_tuned |>
  arrange(MAE) |>
  slice(1)

fig_model_comparison <- ggplot(
  model_comparison_tuned,
  aes(x = reorder(model, MAE), y = MAE, fill = model)
) +
  geom_col(width = 0.7) +
  coord_flip() +
  labs(
    subtitle = "Lower values indicate better prediction accuracy.",
    x = NULL,
    y = "Mean absolute error in 2023 claims"
  ) +
  scale_fill_manual(values = c(
    "Log-linear" = "#4C78A8",
    "Two-part" = "#F58518",
    "Random Forest" = "#54A24B",
    "Random Forest Tuned" = "#2E7D32",
    "XGBoost" = "#E45756",
    "XGBoost Tuned" = "#B279A2"
  )) +
  theme_minimal() +
  theme(legend.position = "none")


# ------------------------------------------------------------
# 10. Child/dependant sensitivity
# ------------------------------------------------------------

df_children <- meps |>
  transmute(
    id = DUPERSID,
    age = as.numeric(AGE1X),
    sex = as.numeric(SEX),
    physical_health = as.numeric(RTHLTH1),
    mental_health = as.numeric(MNHLTH1),
    diabetes = as.numeric(DIABDXY1_M18),
    asthma = as.numeric(ASTHDXY1),
    insurance = as.numeric(INSCOVY1),
    claims_2022 = as.numeric(TOTEXPY1),
    claims_2023 = as.numeric(TOTEXPY2),
    IP_2022 = as.numeric(IPTEXPY1),
    IP_2023 = as.numeric(IPTEXPY2),
    OP_2022 = as.numeric(OPTEXPY1),
    OP_2023 = as.numeric(OPTEXPY2)
  ) |>
  mutate(
    age = ifelse(age < 0, NA, age),
    sex = ifelse(sex < 0, NA, sex),
    physical_health = ifelse(physical_health < 0, NA, physical_health),
    mental_health = ifelse(mental_health < 0, NA, mental_health),
    diabetes = ifelse(diabetes < 0, NA, diabetes),
    asthma = ifelse(asthma < 0, NA, asthma),
    insurance = ifelse(insurance < 0, NA, insurance)
  ) |>
  filter(
    !is.na(age),
    !is.na(claims_2022),
    !is.na(claims_2023),
    claims_2022 >= 0,
    claims_2023 >= 0
  )

child_data_check <- tibble(
  adults = sum(df_children$age >= 18),
  children = sum(df_children$age < 18),
  min_age = min(df_children$age),
  max_age = max(df_children$age)
)

simulate_family_scheme <- function(size, sim_id) {
  adults <- df_children |>
    filter(age >= 18) |>
    slice_sample(n = round(size * 0.8), replace = TRUE)
  
  children <- df_children |>
    filter(age < 18) |>
    slice_sample(n = size - round(size * 0.8), replace = TRUE)
  
  members <- bind_rows(adults, children)
  
  tibble(
    sim_id = sim_id,
    scheme_size = size,
    claims_2022_pm = mean(members$claims_2022),
    claims_2023_pm = mean(members$claims_2023),
    IP_2022_pm = mean(members$IP_2022),
    IP_2023_pm = mean(members$IP_2023),
    OP_2022_pm = mean(members$OP_2022),
    OP_2023_pm = mean(members$OP_2023)
  )
}

family_sim_results <- expand.grid(
  scheme_size = scheme_sizes,
  sim_id = 1:n_sims
) |>
  rowwise() |>
  do(simulate_family_scheme(size = .$scheme_size, sim_id = .$sim_id)) |>
  ungroup()

adult_manual_2022 <- mean(df_children$claims_2022[df_children$age >= 18])
child_manual_2022 <- mean(df_children$claims_2022[df_children$age < 18])
adult_manual_2023 <- mean(df_children$claims_2023[df_children$age >= 18])
child_manual_2023 <- mean(df_children$claims_2023[df_children$age < 18])

family_manual_2022 <- 0.8 * adult_manual_2022 + 0.2 * child_manual_2022
family_manual_2023 <- 0.8 * adult_manual_2023 + 0.2 * child_manual_2023
family_trend_factor <- family_manual_2023 / family_manual_2022

family_pricing_results <- family_sim_results |>
  mutate(
    Z = scheme_size / (scheme_size + k_final_total),
    price_manual = floor_price(family_manual_2022 * family_trend_factor),
    price_experience = floor_price(claims_2022_pm * family_trend_factor),
    price_credibility = floor_price((Z * claims_2022_pm + (1 - Z) * family_manual_2022) * family_trend_factor),
    abs_error_manual = abs(price_manual - claims_2023_pm),
    abs_error_experience = abs(price_experience - claims_2023_pm),
    abs_error_credibility = abs(price_credibility - claims_2023_pm)
  )

family_method_summary <- family_pricing_results |>
  group_by(scheme_size) |>
  summarise(
    Manual_Rating = mean(abs_error_manual),
    Full_Experience = mean(abs_error_experience),
    Credibility = mean(abs_error_credibility),
    .groups = "drop"
  ) |>
  mutate(scenario = "Family mix: 80% adult / 20% child")

adult_method_summary <- method_summary |>
  select(scheme_size, method, mean_abs_error) |>
  pivot_wider(names_from = method, values_from = mean_abs_error) |>
  rename(
    Manual_Rating = `Manual Rating`,
    Full_Experience = `Full Experience`
  ) |>
  mutate(scenario = "Adult-only")

scenario_comparison <- bind_rows(adult_method_summary, family_method_summary) |>
  select(scenario, scheme_size, Manual_Rating, Full_Experience, Credibility)

scenario_long <- scenario_comparison |>
  pivot_longer(
    cols = c(Manual_Rating, Full_Experience, Credibility),
    names_to = "method",
    values_to = "mean_abs_error"
  ) |>
  mutate(
    method = case_when(
      method == "Manual_Rating" ~ "Manual Rating",
      method == "Full_Experience" ~ "Full Experience",
      method == "Credibility" ~ "Credibility"
    )
  )

family_improvement <- family_method_summary |>
  mutate(
    credibility_vs_manual_improvement = (Manual_Rating - Credibility) / Manual_Rating,
    credibility_vs_experience_improvement = (Full_Experience - Credibility) / Full_Experience
  )

fig_family_comparison <- ggplot(scenario_long, aes(x = factor(scheme_size), y = mean_abs_error, fill = method)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = method_colours) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  ) +
  facet_wrap(~ scenario) +
  labs(
    x = "Scheme size",
    y = "Mean absolute pricing error per member",
    fill = "Pricing method"
  )


# ------------------------------------------------------------
# 11. Profit/loss analysis
# ------------------------------------------------------------

profit_results <- pricing_results |>
  mutate(
    premium_manual = price_manual * (1 + target_loading),
    premium_experience = price_experience * (1 + target_loading),
    premium_credibility = price_credibility * (1 + target_loading),
    margin_manual = (premium_manual - claims_2023_pm) / premium_manual,
    margin_experience = (premium_experience - claims_2023_pm) / premium_experience,
    margin_credibility = (premium_credibility - claims_2023_pm) / premium_credibility,
    loss_making_manual = margin_manual < 0,
    loss_making_experience = margin_experience < 0,
    loss_making_credibility = margin_credibility < 0
  )

profit_summary <- profit_results |>
  group_by(scheme_size) |>
  summarise(
    avg_margin_manual = mean(margin_manual),
    avg_margin_experience = mean(margin_experience),
    avg_margin_credibility = mean(margin_credibility),
    prob_loss_manual = mean(loss_making_manual),
    prob_loss_experience = mean(loss_making_experience),
    prob_loss_credibility = mean(loss_making_credibility),
    sd_margin_manual = sd(margin_manual),
    sd_margin_experience = sd(margin_experience),
    sd_margin_credibility = sd(margin_credibility),
    .groups = "drop"
  )

commercial_summary <- profit_summary |>
  select(
    scheme_size,
    prob_loss_manual,
    prob_loss_experience,
    prob_loss_credibility,
    sd_margin_manual,
    sd_margin_experience,
    sd_margin_credibility
  )

loss_prob_long <- commercial_summary |>
  select(
    scheme_size,
    prob_loss_manual,
    prob_loss_experience,
    prob_loss_credibility
  ) |>
  pivot_longer(
    cols = starts_with("prob_loss"),
    names_to = "method",
    values_to = "prob_loss"
  ) |>
  mutate(
    method = case_when(
      method == "prob_loss_manual" ~ "Manual Rating",
      method == "prob_loss_experience" ~ "Full Experience",
      method == "prob_loss_credibility" ~ "Credibility"
    )
  )

fig_loss_probability <- ggplot(
  loss_prob_long,
  aes(x = factor(scheme_size), y = prob_loss, fill = method)
) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = method_colours) +
  scale_y_continuous(labels = scales::percent) +
  labs(
    subtitle = "Lower values indicate lower risk of underpricing.",
    x = "Scheme size",
    y = "Probability of loss",
    fill = "Pricing method"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )
# ------------------------------------------------------------
# 12. Stress testing
# ------------------------------------------------------------

inflation_stress <- pricing_results |>
  mutate(
    claims_2023_stress = claims_2023_pm * 1.10,
    abs_error_manual_stress = abs(price_manual - claims_2023_stress),
    abs_error_experience_stress = abs(price_experience - claims_2023_stress),
    abs_error_credibility_stress = abs(price_credibility - claims_2023_stress)
  ) |>
  group_by(scheme_size) |>
  summarise(
    manual_MAE = mean(abs_error_manual_stress),
    experience_MAE = mean(abs_error_experience_stress),
    credibility_MAE = mean(abs_error_credibility_stress),
    .groups = "drop"
  )

large_ip_claim <- quantile(df_clean$IP_2023[df_clean$IP_2023 > 0], 0.99, na.rm = TRUE)

ip_shock <- pricing_results |>
  mutate(
    claims_2023_stress = claims_2023_pm + large_ip_claim / scheme_size,
    abs_error_manual_stress = abs(price_manual - claims_2023_stress),
    abs_error_experience_stress = abs(price_experience - claims_2023_stress),
    abs_error_credibility_stress = abs(price_credibility - claims_2023_stress)
  ) |>
  group_by(scheme_size) |>
  summarise(
    manual_MAE = mean(abs_error_manual_stress),
    experience_MAE = mean(abs_error_experience_stress),
    credibility_MAE = mean(abs_error_credibility_stress),
    .groups = "drop"
  )


# ------------------------------------------------------------
# 13. Save outputs
# ------------------------------------------------------------

final_improvement <- method_summary |>
  select(scheme_size, method, mean_abs_error) |>
  pivot_wider(names_from = method, values_from = mean_abs_error) |>
  mutate(
    credibility_vs_manual_improvement = (`Manual Rating` - Credibility) / `Manual Rating`,
    credibility_vs_experience_improvement = (`Full Experience` - Credibility) / `Full Experience`
  )


# Tables

write.csv(claims_summary, "outputs/tables/01_claims_summary.csv", row.names = FALSE)
write.csv(mixed_profiles, "outputs/tables/02_mixed_profiles.csv", row.names = FALSE)
write.csv(cred_sensitivity_overall, "outputs/tables/03_credibility_sensitivity_overall.csv", row.names = FALSE)
write.csv(method_summary, "outputs/tables/04_method_summary.csv", row.names = FALSE)
write.csv(final_improvement, "outputs/tables/05_final_improvement.csv", row.names = FALSE)
write.csv(mixed_method_summary, "outputs/tables/06_mixed_method_summary.csv", row.names = FALSE)
write.csv(split_summary, "outputs/tables/07_ip_op_summary.csv", row.names = FALSE)
write.csv(best_split_k, "outputs/tables/08_best_ip_op_k.csv", row.names = FALSE)
write.csv(model_comparison_tuned, "outputs/tables/09_model_comparison_tuned.csv", row.names = FALSE)
write.csv(rf_importance, "outputs/tables/10_rf_importance.csv", row.names = FALSE)
write.csv(xgb_importance, "outputs/tables/11_xgb_importance.csv", row.names = FALSE)
write.csv(scenario_comparison, "outputs/tables/12_child_sensitivity_comparison.csv", row.names = FALSE)
write.csv(profit_summary, "outputs/tables/13_profit_summary.csv", row.names = FALSE)
write.csv(inflation_stress, "outputs/tables/14_inflation_stress.csv", row.names = FALSE)
write.csv(ip_shock, "outputs/tables/15_ip_shock.csv", row.names = FALSE)


# Figures
ggsave("outputs/figures/01_claims_distribution.png", fig_claims_distribution, width = 8, height = 5, dpi = 300)
ggsave("outputs/figures/02_scheme_volatility.png", fig_scheme_volatility, width = 7, height = 5, dpi = 300)
ggsave("outputs/figures/03_scheme_claim_range.png",fig_scheme_claim_range, width = 7, height = 5, dpi = 300)
ggsave("outputs/figures/04_pricing_accuracy.png", fig_pricing_accuracy, width = 8, height = 5, dpi = 300)
ggsave("outputs/figures/05_mixed_pricing_accuracy.png", fig_mixed_pricing_accuracy, width = 9, height = 5, dpi = 300)
ggsave("outputs/figures/06_ip_op_accuracy.png", fig_split_accuracy, width = 9, height = 5, dpi = 300)
ggsave("outputs/figures/07_ip_op_k_sensitivity.png", fig_split_k, width = 11, height = 7, dpi = 300)
ggsave("outputs/figures/08_model_comparison.png", fig_model_comparison, width = 8, height = 5, dpi = 300)
ggsave("outputs/figures/09_child_sensitivity.png", fig_family_comparison, width = 9, height = 5, dpi = 300)
ggsave("outputs/figures/10_loss_probability.png", fig_loss_probability, width = 8, height = 5, dpi = 300)

save(
  claims_summary,
  cred_sensitivity_overall,
  mixed_profiles,
  method_summary,
  final_improvement,
  mixed_method_summary,
  split_summary,
  model_comparison_tuned,
  rf_importance,
  xgb_importance,
  scenario_comparison,
  profit_summary,
  inflation_stress,
  ip_shock,
  
  fig_claims_distribution,
  fig_scheme_volatility,
  fig_scheme_claim_range,
  fig_pricing_accuracy,
  fig_mixed_pricing_accuracy,
  fig_split_accuracy,
  fig_split_k,
  fig_model_comparison,
  fig_family_comparison,
  fig_loss_probability,
  
  file = "thesis_outputs.RData"
)
