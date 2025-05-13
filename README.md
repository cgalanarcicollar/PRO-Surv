# PRO-Surv

Code and materials for the paper:  
[**Patient-reported outcomes and survival analysis of chronic obstructive pulmonary disease patients: A two-stage joint modelling approach**](https://raco.cat/index.php/SORT/article/view/431617)  
Galán-Arcicollar, C., Najera-Zuloaga, J., & Lee, D.-J. (2024). *SORT*, 48(2), 155–182. https://doi.org/10.57645/20.8080.02.17

## Overview

This repository contains the code used to perform the statistical analysis described in our paper. We propose a two-stage joint modeling including beta-binomial distribution to analyze the relationship between patient-reported outcomes (PROs) and survival outcomes.
It contains the R code for all relevant analyses in the paper. Real data and settings from the first scenario simulation (based on real data), are excluded because the data cannot be made publicly available.

##  Main Features

- **Beta-Binomial longitudinal modeling**: Incorporates a flexible and realistic modeling approach for patient-reported outcomes (PROs), which are often bounded and overdispersed.
  
- **Two-stage joint modeling framework**: Uses a computationally efficient frequentist two-stage method. The longitudinal model is fitted first, and its predicted values are then included as covariates in the survival model.
  
- **Comparison of modeling approaches**:
  - The proposed **Beta-Binomial + two-stage joint model**.
  - A **classical joint model** assuming a Gaussian longitudinal component.
  - An **extended (time-varying) Cox model** that includes PROs directly as time-dependent covariates.
  
- **Simulation study**: Assesses and compares the performance and bias of the different approaches under various scenarios.
