
# **📊 DEMO Design: Dose Exploration, Monitoring, and Optimization using Biological and Clinical Outcomes**

<!-- badges: start -->
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R Version](https://img.shields.io/badge/R->=4.2-blue)](https://cran.r-project.org/)
[![JAGS Required](https://img.shields.io/badge/JAGS-Required-red)](http://mcmc-jags.sourceforge.net/)
<!-- badges: end -->

## Overview

**DEMOdesign** is an R package implementing the **Dose Exploration, Monitoring, and Optimization (DEMO)** design—a Bayesian adaptive Phase I/II trial framework tailored for oncology dose-finding, particularly for immunotherapies and targeted therapies. It selects an **Optimal Therapeutic Dose (OTD)** that maximizes **Restricted Mean Survival Time (RMST)** while ensuring safety, biological activity, and clinical efficacy across three stages:

1. **Exploration**: Identifies safe, biologically active doses using biomarkers ($Y_B$) and toxicity ($Y_T$).
2. **Monitoring**: Refines doses based on response ($Y_R$) and additional safety data.
3. **Optimization**: Selects the OTD using long-term survival ($Y_S$) via a Weibull model.

This package supports:

> Yang, C.-H., Thall, P. F., & Lin, R. (2024). **DEMO: Dose Exploration, Monitoring, and Optimization Using Biological and Clinical Outcomes**. *Annals of Applied Statistics*. (Under Review).

---

## Features

- **Three-Stage Bayesian Design**:
  - **Stage 1**: BOIN-based exploration with biomarker and toxicity screening.
  - **Stage 2**: Monitors admissible doses for efficacy and safety.
  - **Stage 3**: Optimizes dose selection with RMST.
- **Endpoints**: Integrates `Y_B` (biomarker), `Y_T` (toxicity), `Y_R` (response), `Y_S` (survival).

---

## Installation

Install from GitHub:

``` r
remotes::install_github("cyang728/DEMOdesign")
```

## Example

The following example demonstrates how to use DEMOdesign to simulate a Bayesian adaptive dose-finding trial.

``` r
library(DEMOdesign)

# Run DEMO_design with example data
result <- DEMO_design(seed=1,
                      doses = c(0.05, 0.10, 0.20, 0.45, 0.65, 0.85),
                      Y_B_sim = c(2.00, 2.01, 2.08, 2.76, 3.75, 4.73),  # biomarker means
                      sigma2_B_sim = 1, 
                      Y_T_sim = c(0.01, 0.02, 0.03, 0.06, 0.13, 0.26),  # toxicity rates
                      Y_R_sim = c(0.04, 0.05, 0.08, 0.20, 0.35, 0.47),  # response rates
                      lambdaT_sim = c(0.8, 0.6, 0.6, 0.25, 0.2, 0.1),   # scale for Weibull survival
                      shape_sim = 1.5,      # Weibull shape parameter
                      delta1_sim = 3,       # effect of toxicity on survival
                      delta2_sim = -2,      # effect of response on survival
                      censored_time = 24    # administrative censoring time 
                      )

# View the optimal therapeutic dose (OTD)
print(result$OTD)

# View summary of patient-level trial data
head(result$trial)
```

## Function Documentation

#### Inputs

- `seed`: Integer. Random seed for reproducibility.
- `doses`: Numeric vector. Dose levels to evaluate.
- `Y_B_sim`: Numeric vector. Mean biomarker values for each dose.
- `sigma2_B_sim`: Numeric. Variance of the biomarker outcome.
- `Y_T_sim`: Numeric vector. Toxicity probabilities for each dose.
- `Y_R_sim`: Numeric vector. Response probabilities for each dose.
- `lambdaT_sim`: Numeric vector. Weibull scale parameters for survival at each dose.
- `shape_sim`: Numeric. Weibull shape parameter.
- `delta1_sim`: Numeric. Effect of toxicity on survival time.
- `delta2_sim`: Numeric. Effect of response on survival time.
- `censored_time`: Numeric. Administrative censoring time (e.g., 24 months).

#### Returns

A list with the following elements:

- `OTD`: The selected Optimal Therapeutic Dose.
- `trial`: A data frame with patient-level trial data, including assigned dose, biomarker, toxicity, response, and survival outcomes.
- `summary`: (Optional) Summary statistics across the trial (if implemented).

