# ============================================================
# MSc Actuarial Science Thesis
# Covid-period robustness check
# Credibility-Based Pricing of Corporate Health Insurance Schemes
# Data: MEPS HC-245, rolling one-year tests
# ============================================================

# ------------------------------------------------------------
# 0. Packages and settings
# ------------------------------------------------------------

library(haven)
library(dplyr)
library(tidyr)
library(ggplot2)

method_colours <- c(
  "Credibility" = "#619CFF",
  "Manual Rating" = "#F8766D",
  "Full Experience" = "#00BA38"
)

# Set seed so simulations can be reproduced
set.seed(183)

scheme_sizes <- c(10, 50, 100, 500, 1500)
age_profiles <- c("Young", "Middle", "Older")
risk_profiles <- c("Low", "Medium", "High")
n_sims <- 500
k_values <- seq(5, 75, by = 5)

floor_price <- function(x) pmax(x, 1)

# Output folders
if (!dir.exists("outputs")) dir.create("outputs")
if (!dir.exists("outputs/covid")) dir.create("outputs/covid")
if (!dir.exists("outputs/covid/tables")) dir.create("outputs/covid/tables", recursive = TRUE)
if (!dir.exists("outputs/covid/figures")) dir.create("outputs/covid/figures", recursive = TRUE)


# ------------------------------------------------------------
# 1. Load data
# ------------------------------------------------------------

meps245 <- read_dta("h245.dta")

raw_data_check <- tibble(
  file = "HC-245",
  rows = nrow(meps245),
  columns = ncol(meps245)
)


# ------------------------------------------------------------
# 2. Define HC-245 rolling year pairs
# ------------------------------------------------------------

year_pairs <- tibble(
  period = c("2019 to 2020", "2020 to 2021", "2021 to 2022"),
  start_year = c(2019, 2020, 2021),
  end_year = c(2020, 2021, 2022),
  start_num = c(1, 2, 3),
  end_num = c(2, 3, 4)
)


# ------------------------------------------------------------
# 3. Build clean adult dataset for one period
# ------------------------------------------------------------

build_period_data <- function(data, start_num, end_num, period_label) {
  
  age_var <- paste0("AGEY", start_num, "X")
  physical_var <- paste0("RTHLTH", start_num)
  mental_var <- paste0("MNHLTH", start_num)
  diabetes_var <- paste0("DIABDXY", start_num, "_M18")
  asthma_var <- paste0("ASTHDXY", start_num)
  insurance_var <- paste0("INSCOVY", start_num)
  
  claims_start_var <- paste0("TOTEXPY", start_num)
  claims_end_var <- paste0("TOTEXPY", end_num)
  
  required_vars <- c(
    "DUPERSID", "SEX",
    age_var, physical_var, mental_var, diabetes_var, asthma_var, insurance_var,
    claims_start_var, claims_end_var
  )
  
  missing_vars <- setdiff(required_vars, names(data))
  
  if (length(missing_vars) > 0) {
    stop(paste("Missing variables:", paste(missing_vars, collapse = ", ")))
  }
  
  tibble(
    period = period_label,
    id = as.character(data$DUPERSID),
    age = as.numeric(data[[age_var]]),
    sex = as.numeric(data$SEX),
    physical_health = as.numeric(data[[physical_var]]),
    mental_health = as.numeric(data[[mental_var]]),
    diabetes = as.numeric(data[[diabetes_var]]),
    asthma = as.numeric(data[[asthma_var]]),
    insurance = as.numeric(data[[insurance_var]]),
    claims_start = as.numeric(data[[claims_start_var]]),
    claims_end = as.numeric(data[[claims_end_var]])
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
      !is.na(claims_start),
      !is.na(claims_end),
      claims_start >= 0,
      claims_end >= 0
    )
}


# ------------------------------------------------------------
# 4. Run credibility analysis for one period
# ------------------------------------------------------------

run_period_analysis <- function(df, period_label) {
  
  clean_data_check <- tibble(
    period = period_label,
    rows = nrow(df),
    columns = ncol(df),
    min_age = min(df$age, na.rm = TRUE),
    max_age = max(df$age, na.rm = TRUE),
    mean_claims_start = mean(df$claims_start, na.rm = TRUE),
    mean_claims_end = mean(df$claims_end, na.rm = TRUE)
  )
  
  claims_summary <- tibble(
    period = period_label,
    year_position = c("Start year", "End year"),
    mean_claims = c(mean(df$claims_start, na.rm = TRUE), mean(df$claims_end, na.rm = TRUE)),
    median_claims = c(median(df$claims_start, na.rm = TRUE), median(df$claims_end, na.rm = TRUE)),
    p95_claims = c(
      quantile(df$claims_start, 0.95, na.rm = TRUE),
      quantile(df$claims_end, 0.95, na.rm = TRUE)
    ),
    max_claims = c(max(df$claims_start, na.rm = TRUE), max(df$claims_end, na.rm = TRUE))
  )
  
  # Risk score model
  df_model <- df |>
    mutate(
      claims_start_log = log(claims_start + 1),
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
      !is.na(claims_start_log)
    ) |>
    droplevels()
  
  candidate_predictors <- c(
    "age", "sex", "physical_health", "mental_health",
    "diabetes", "asthma", "insurance"
  )
  
  usable_predictors <- candidate_predictors[
    sapply(df_model[candidate_predictors], function(x) {
      if (is.factor(x)) {
        nlevels(droplevels(x)) >= 2
      } else {
        length(unique(na.omit(x))) >= 2
      }
    })
  ]
  
  risk_formula <- as.formula(
    paste("claims_start_log ~", paste(usable_predictors, collapse = " + "))
  )
  
  risk_model <- lm(risk_formula, data = df_model)
  
  df_clean <- df_model |>
    mutate(
      risk_score_model = exp(predict(risk_model)) - 1,
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
    group_by(period, age_profile, risk_profile_model) |>
    summarise(
      members = n(),
      avg_risk_score = mean(risk_score_model, na.rm = TRUE),
      avg_claims_start = mean(claims_start, na.rm = TRUE),
      avg_claims_end = mean(claims_end, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Scheme simulation
  simulate_scheme <- function(size, age_prof, risk_prof, sim_id) {
    pool <- df_clean |>
      filter(
        age_profile == age_prof,
        risk_profile_model == risk_prof
      )
    
    if (nrow(pool) == 0) {
      stop(paste("No members found for", period_label, age_prof, risk_prof))
    }
    
    sample_members <- pool |>
      slice_sample(n = size, replace = TRUE)
    
    tibble(
      period = period_label,
      sim_id = sim_id,
      scheme_size = size,
      age_profile = age_prof,
      risk_profile = risk_prof,
      claims_start_pm = mean(sample_members$claims_start, na.rm = TRUE),
      claims_end_pm = mean(sample_members$claims_end, na.rm = TRUE)
    )
  }
  
  sim_grid <- expand.grid(
    scheme_size = scheme_sizes,
    age_profile = age_profiles,
    risk_profile = risk_profiles,
    sim_id = 1:n_sims,
    stringsAsFactors = FALSE
  )
  
  sim_results <- bind_rows(lapply(seq_len(nrow(sim_grid)), function(i) {
    simulate_scheme(
      size = sim_grid$scheme_size[i],
      age_prof = sim_grid$age_profile[i],
      risk_prof = sim_grid$risk_profile[i],
      sim_id = sim_grid$sim_id[i]
    )
  }))
  
  scheme_summary <- sim_results |>
    group_by(period, scheme_size, age_profile, risk_profile) |>
    summarise(
      avg_start_claims = mean(claims_start_pm, na.rm = TRUE),
      avg_end_claims = mean(claims_end_pm, na.rm = TRUE),
      sd_end_claims = sd(claims_end_pm, na.rm = TRUE),
      p10_end_claims = quantile(claims_end_pm, 0.10, na.rm = TRUE),
      median_end_claims = median(claims_end_pm, na.rm = TRUE),
      p90_end_claims = quantile(claims_end_pm, 0.90, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Pricing methods
  market_start <- mean(df_clean$claims_start, na.rm = TRUE)
  market_end <- mean(df_clean$claims_end, na.rm = TRUE)
  trend_factor <- market_end / market_start
  
  manual_rates <- df_clean |>
    group_by(age_profile, risk_profile_model) |>
    summarise(
      manual_start = mean(claims_start, na.rm = TRUE),
      .groups = "drop"
    ) |>
    rename(risk_profile = risk_profile_model)
  
  apply_total_pricing <- function(k_value) {
    sim_results |>
      left_join(manual_rates, by = c("age_profile", "risk_profile")) |>
      mutate(
        k = k_value,
        Z = scheme_size / (scheme_size + k),
        price_manual = floor_price(manual_start * trend_factor),
        price_experience = floor_price(claims_start_pm * trend_factor),
        price_credibility = floor_price((Z * claims_start_pm + (1 - Z) * manual_start) * trend_factor),
        error_manual = price_manual - claims_end_pm,
        error_experience = price_experience - claims_end_pm,
        error_credibility = price_credibility - claims_end_pm,
        abs_error_manual = abs(error_manual),
        abs_error_experience = abs(error_experience),
        abs_error_credibility = abs(error_credibility),
        lr_manual = claims_end_pm / price_manual,
        lr_experience = claims_end_pm / price_experience,
        lr_credibility = claims_end_pm / price_credibility
      )
  }
  
  make_pricing_long <- function(pricing_data) {
    pricing_data |>
      select(
        period, scheme_size, age_profile, risk_profile, sim_id,
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
  
  evaluate_credibility_k <- function(k_value) {
    apply_total_pricing(k_value) |>
      group_by(period, scheme_size) |>
      summarise(
        k = k_value,
        mean_abs_error_credibility = mean(abs_error_credibility, na.rm = TRUE),
        median_abs_error_credibility = median(abs_error_credibility, na.rm = TRUE),
        .groups = "drop"
      )
  }
  
  cred_sensitivity <- bind_rows(lapply(k_values, evaluate_credibility_k))
  
  cred_sensitivity_overall <- cred_sensitivity |>
    group_by(period, k) |>
    summarise(
      overall_mean_abs_error = mean(mean_abs_error_credibility, na.rm = TRUE),
      overall_median_abs_error = mean(median_abs_error_credibility, na.rm = TRUE),
      .groups = "drop"
    )
  
  best_total_k <- cred_sensitivity_overall |>
    arrange(overall_mean_abs_error) |>
    slice(1)
  
  k_final_total <- best_total_k$k
  
  pricing_results <- apply_total_pricing(k_final_total)
  pricing_long <- make_pricing_long(pricing_results)
  
  method_summary <- pricing_long |>
    group_by(period, scheme_size, method) |>
    summarise(
      mean_abs_error = mean(abs_error, na.rm = TRUE),
      median_abs_error = median(abs_error, na.rm = TRUE),
      sd_abs_error = sd(abs_error, na.rm = TRUE),
      p5_abs_error = quantile(abs_error, 0.05, na.rm = TRUE),
      p95_abs_error = quantile(abs_error, 0.95, na.rm = TRUE),
      .groups = "drop"
    )
  
  loss_ratio_summary <- pricing_results |>
    group_by(period, scheme_size, age_profile, risk_profile) |>
    summarise(
      avg_LR_manual = mean(lr_manual, na.rm = TRUE),
      avg_LR_experience = mean(lr_experience, na.rm = TRUE),
      avg_LR_credibility = mean(lr_credibility, na.rm = TRUE),
      .groups = "drop"
    )
 
  final_improvement <- method_summary |>
    select(period, scheme_size, method, mean_abs_error) |>
    pivot_wider(names_from = method, values_from = mean_abs_error) |>
    mutate(
      credibility_vs_manual_improvement = (`Manual Rating` - Credibility) / `Manual Rating`,
      credibility_vs_experience_improvement = (`Full Experience` - Credibility) / `Full Experience`
    )
  
  list(
    clean_data_check = clean_data_check,
    claims_summary = claims_summary,
    risk_band_check = risk_band_check,
    scheme_summary = scheme_summary,
    cred_sensitivity = cred_sensitivity,
    cred_sensitivity_overall = cred_sensitivity_overall,
    best_total_k = best_total_k,
    method_summary = method_summary,
    loss_ratio_summary = loss_ratio_summary,
    final_improvement = final_improvement
  )
}


# ------------------------------------------------------------
# 5. Run rolling Covid-period tests
# ------------------------------------------------------------

covid_runs <- lapply(seq_len(nrow(year_pairs)), function(i) {
  period_data <- build_period_data(
    data = meps245,
    start_num = year_pairs$start_num[i],
    end_num = year_pairs$end_num[i],
    period_label = year_pairs$period[i]
  )
  
  run_period_analysis(
    df = period_data,
    period_label = year_pairs$period[i]
  )
})


# ------------------------------------------------------------
# 6. Combine outputs
# ------------------------------------------------------------

covid_clean_data_check <- bind_rows(lapply(covid_runs, function(x) x$clean_data_check))
covid_claims_summary <- bind_rows(lapply(covid_runs, function(x) x$claims_summary))
covid_risk_band_check <- bind_rows(lapply(covid_runs, function(x) x$risk_band_check))
covid_scheme_summary <- bind_rows(lapply(covid_runs, function(x) x$scheme_summary))
covid_cred_sensitivity <- bind_rows(lapply(covid_runs, function(x) x$cred_sensitivity))
covid_cred_sensitivity_overall <- bind_rows(lapply(covid_runs, function(x) x$cred_sensitivity_overall))
covid_best_k <- bind_rows(lapply(covid_runs, function(x) x$best_total_k))
covid_method_summary <- bind_rows(lapply(covid_runs, function(x) x$method_summary))
covid_loss_ratio_summary <- bind_rows(lapply(covid_runs, function(x) x$loss_ratio_summary))
covid_final_improvement <- bind_rows(lapply(covid_runs, function(x) x$final_improvement))

covid_best_methods <- covid_method_summary |>
  group_by(period, scheme_size) |>
  slice_min(mean_abs_error, n = 1, with_ties = FALSE) |>
  ungroup() |>
  select(period, scheme_size, method, mean_abs_error)

covid_period_summary <- covid_best_methods |>
  group_by(period) |>
  summarise(
    best_method_size_10 = method[scheme_size == 10][1],
    best_method_size_50 = method[scheme_size == 50][1],
    best_method_size_100 = method[scheme_size == 100][1],
    best_method_size_500 = method[scheme_size == 500][1],
    best_method_size_1500 = method[scheme_size == 1500][1],
    .groups = "drop"
  )


# ------------------------------------------------------------
# 7. Figures
# ------------------------------------------------------------

fig_covid_pricing_accuracy <- ggplot(
  covid_method_summary,
  aes(x = factor(scheme_size), y = mean_abs_error, fill = method)
) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = method_colours) +
  facet_wrap(~ period, nrow = 1) +
  scale_y_continuous(labels = scales::comma) +
  labs(
    subtitle = "Rolling one-year tests using MEPS HC-245. Lower values indicate better pricing accuracy.",
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

fig_covid_k_sensitivity <- ggplot(
  covid_cred_sensitivity_overall,
  aes(x = k, y = overall_mean_abs_error)
) +
  geom_line() +
  geom_point() +
  facet_wrap(~ period, scales = "free_y") +
  labs(
    title = "Covid-Period Sensitivity to Credibility Parameter",
    subtitle = "Lower values indicate better credibility pricing accuracy.",
    x = "Credibility parameter k",
    y = "Overall mean absolute pricing error"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold")
  )


# ------------------------------------------------------------
# 8. Save outputs
# ------------------------------------------------------------

write.csv(covid_claims_summary, "outputs/covid/tables/01_covid_claims_summary.csv", row.names = FALSE)
write.csv(covid_best_k, "outputs/covid/tables/02_covid_best_k.csv", row.names = FALSE)
write.csv(covid_method_summary, "outputs/covid/tables/03_covid_method_summary.csv", row.names = FALSE)
write.csv(covid_final_improvement, "outputs/covid/tables/04_covid_final_improvement.csv", row.names = FALSE)
write.csv(covid_best_methods, "outputs/covid/tables/05_covid_best_methods.csv", row.names = FALSE)
write.csv(covid_period_summary, "outputs/covid/tables/06_covid_period_summary.csv", row.names = FALSE)
write.csv(covid_cred_sensitivity_overall, "outputs/covid/tables/07_covid_credibility_sensitivity_overall.csv", row.names = FALSE)

ggsave("outputs/covid/figures/01_covid_pricing_accuracy.png", fig_covid_pricing_accuracy, width = 11, height = 5, dpi = 300)

save(
  covid_best_k,
  covid_method_summary,
  covid_best_methods,
  covid_period_summary,
  covid_cred_sensitivity_overall,
  fig_covid_pricing_accuracy,
  file = "covid_outputs.RData"
)
