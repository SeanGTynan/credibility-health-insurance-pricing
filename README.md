# Credibility-Based Pricing of Corporate Health Insurance Schemes

This repository contains the R code used for my MSc Actuarial Science thesis.

The project examines whether credibility pricing can improve the accuracy and stability of corporate health insurance scheme pricing compared with manual rating and full experience rating.

## Repository Structure

- `code/Main_Analysis.R` – main pricing analysis and robustness checks
- `code/Covid_Robustness.R` – COVID-19 period robustness analysis
- `data/README.md` – information on the MEPS datasets used

## Data

The analysis uses publicly available Medical Expenditure Panel Survey (MEPS) data.

- HC-252 Panel 27 Longitudinal File – main 2022–2023 analysis
- HC-245 Panel 24 Longitudinal File – COVID-19 period robustness analysis

The raw MEPS data files are not included in this repository.

## Software

The analysis was carried out in R using the following packages:

- haven
- dplyr
- tidyr
- ggplot2
- ranger
- xgboost
