
# DEMOdesign: Bayesian Adaptive Dose-Finding Using Biological and Clinical Outcomes

<!-- badges: start -->
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R Version](https://img.shields.io/badge/R->=4.2-blue)](https://cran.r-project.org/)
[![JAGS Required](https://img.shields.io/badge/JAGS-Required-red)](http://mcmc-jags.sourceforge.net/)
<!-- badges: end -->

## Overview

**DEMOdesign** implements the **Dose Exploration, Monitoring, and Optimization (DEMO) design**, a **Bayesian adaptive phase I/II trial design** that integrates:
- **Early biological responses (YB)** as a biomarker for dose activity.
- **Short-term toxicity (YT)** for safety monitoring.
- **Intermediate response (YR)** for clinical efficacy assessment.
- **Long-term survival (YS)** modeled via a **Weibull survival model**.

The DEMO design optimally selects the dose that maximizes restricted mean survival time (RMST) while ensuring acceptable toxicity, sufficient biological activity, and clinical efficacy.

This package is based on the paper:

> Yang, C.-H., Thall, P. F., & Lin, R. (2024). **DEMO: Dose Exploration, Monitoring, and Optimization Using Biological and Clinical Outcomes**. *Annals of Applied Statistics*. (Under Review).


## **📊 DEMO Design: Three-Stage Bayesian Adaptive Dose-Finding**

The **DEMOdesign** follows a **three-stage Bayesian adaptive framework**, incorporating **biomarker response, short-term toxicity, intermediate response, and survival outcomes** to guide dose selection.

---

### **1️⃣ Stage 1: Dose Exploration (Early Stopping for Safety & Activity)**

📌 **Objective:** Identify **biologically active** and **non-toxic** doses.  
📌 **Method:** Uses **Bayesian Optimal Interval (BOIN)** for **dose escalation/de-escalation**.  
📌 **Decision Rules:**
- **Eliminate doses** that show **high toxicity**.
- **Drop doses** with **insufficient biomarker activity (YB)**.

✔ **Interim Analyses in Stage 1:**  
- **Midpoint Analysis**: Conducted **after 50% of cohorts are enrolled**.
- **Final Exploration Analysis**: Determines biologically active doses at the **end of dose escalation**.

---

### **2️⃣ Stage 2: Dose Monitoring (Randomized Screening for Clinical Efficacy)**

📌 **Objective:** Among biologically active doses, **identify those with acceptable toxicity and sufficient clinical efficacy**.  
📌 **Method:** **Adaptive randomization** among remaining doses to estimate:
- **Short-term toxicity (YT)** and **tumor response (Y_R)**.
- **Early signs of unacceptable efficacy**.

✔ **Interim Analyses in Stage 2:**  
- **Multiple analyses** performed every **3-6 patients per dose**.
- **Stopping rules:**  
  - **Drop doses** with a **high probability of unacceptable toxicity**.
  - **Stop randomization** to doses with **poor response rates**.

---

### **3️⃣ Stage 3: Dose Optimization (Final Selection Using RMST)**

📌 **Objective:** Select the **Optimal Therapeutic Dose (OTD)** that **maximizes long-term survival (YS)**.  
📌 **Method:** Uses a **Bayesian Weibull survival model** to evaluate **Restricted Mean Survival Time (RMST)**.  
📌 **Final Decision Rule:**  

\[
\text{OTD} = \underset{d}{\arg\max} \, \text{RMST}(d) \quad \text{subject to} \quad P(Toxicity_d > c_T) < c_T, \quad P(Response_d < c_R) < c_R
\]

✔ **Interim Analyses in Stage 3:**  
- **First survival analysis** after **50% of patients are enrolled**.  
- **Final selection** of **OTD** at the **end of the study**.

---

## **🔄 Summary: Number of Interim Analyses**
| **Stage**  | **Objective**                                      | **Number of Interim Analyses**  |
|------------|---------------------------------------------------|--------------------------------|
| **Stage 1** | Identify biologically active & safe doses       | **1-2** (Midpoint + Final BOIN) |
| **Stage 2** | Screen doses for efficacy & toxicity            | **Multiple** (Every few cohorts) |
| **Stage 3** | Select OTD based on survival (RMST)             | **2** (Midpoint + Final selection) |

✔ **Total Number of Interim Analyses: ~5-7**  
✔ **Dynamically adjusts dose selection throughout the trial** to ensure the best patient outcomes.

---

This **three-stage Bayesian adaptive framework** ensures that **safe, biologically active, and clinically effective doses** are selected while optimizing **patient survival outcomes**.


## Installation

You can install the development version of DEMOdesign like so:

``` r
remotes::install_github("cyang728/DEMOdesign")
```

## Example

The following example demonstrates how to use DEMOdesign to simulate a Bayesian adaptive dose-finding trial.

``` r
library(DEMOdesign)

# Run DEMO_design with default parameters
result <- DEMO_design(seed=1,
                      Y_B_sim = c(2.00, 2.01, 2.08, 2.76, 3.75, 4.73), 
                      sigma2_B_sim = 1, 
                      Y_T_sim = c(0.01, 0.02, 0.03, 0.06, 0.13, 0.26), 
                      Y_R_sim = c(0.04, 0.05, 0.08, 0.20, 0.35, 0.47),
                      lambdaT_sim = c(0.8, 0.6, 0.6, 0.25, 0.2, 0.1),
                      shape_sim = 1.5,
                      delta1_sim = 3, delta2_sim = -2,
                      censored_time = 24)

# Print the optimal therapeutic dose (OTD)
print(result$OTD)

# View trial summary
head(result$trial)
```




