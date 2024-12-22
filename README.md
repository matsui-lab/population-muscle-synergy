# Unbiased estimation of the population-level motor module

This repository contains a comprehensive pipeline for generating and analyzing muscle synergy data using Non-negative Matrix Factorization (NMF). The pipeline is designed for large-scale simulation and synergy analysis, including individual and population-level computations, synergy matching, and result visualization.

---

## Table of Contents

1. [Overview](#overview)
2. [Features](#features)
3. [Requirements](#requirements)
4. [Usage Instructions](#usage-instructions)
5. [File Descriptions](#file-descriptions)
6. [Visualization Outputs](#visualization-outputs)

---

## Overview


**Background**: Motor modules, or muscle synergies, are fundamental neurophysiological commands essential for muscle coordination. In clinical settings, characterizing these modules at the population level is crucial for understanding pathophysiological mechanisms and assessing intervention outcomes. Traditional approaches estimate individual motor modules and align them across subjects. However, issues such as dependency on initial values, alignment ambiguities, and numerical instability in concatenated non-negative matrix factorization (cNMF) complicate the reliable characterization of the latent population distribution.

**Methods**: To address these challenges, we developed a robust statistical framework aimed at accurately capturing population-level motor module characteristics. We conducted simulation experiments to evaluate the performance of existing methods, including Module Matching and cNMF, under varying sample sizes and noise conditions. Subsequently, we proposed a novel estimation procedure leveraging functional data analysis to estimate motor modules directly from the population mean of muscle activity.

**Results**: Our simulation results indicated that existing methods exhibit significant estimation biases and fail to converge as sample sizes increase, particularly in noisy environments. In contrast, our proposed approach demonstrated superior accuracy and robustness, reducing estimation bias and achieving asymptotic convergence with as few as 25 samples. Additionally, it maintained reliability even with 10% contaminated data, outperforming conventional methods in robustness to outliers.

**Conclusions**: These findings underscore critical limitations in widely used NMF-based approaches and introduce a robust alternative for unbiased estimation of population-level motor modules. Our framework holds significant potential for clinical and experimental applications, facilitating precise evaluation of motor control mechanisms and intervention effects across diverse populations. Future research should focus on refining this approach to accommodate heterogeneous subpopulations and expanding its applicability to various motor tasks.


This project implements a pipeline for:
- Generating synthetic muscle activity data.
- Calculating synergies using NMF for both individual and population levels.
- Matching synergies across individuals.
- Visualizing results for both `W` (activation patterns) and `H` (muscle contributions) matrices.
- Comparing reconstruction errors across different methods (individual, population, concatenation).

The pipeline is fully customizable for various parameters such as the number of synergies, sample sizes, and noise levels.

---

## Features

- **Simulation of Muscle Activity**:
  Generates synthetic muscle activity data with realistic noise and variability.
  
- **Synergy Computation**:
  Uses NMF to calculate individual and population-level synergies.

- **Synergy Matching**:
  Aligns synergies across different samples for comparison.

- **Visualization**:
  Produces publication-ready plots for synergy comparisons in both `W` and `H` matrices.

- **Residual Analysis**:
  Evaluates reconstruction errors using Root Mean Square Error (RMSE) metrics.

---

## Requirements

- **R Version**: >= 4.0
- **R Packages**:
  - `parallel`
  - `data.table`
  - `NMF`
  - `tidyverse`
  - `ggplot2`
  - `gridExtra`
  - `splines`

Ensure all dependencies are installed before running the RMarkdown files.

---

## Usage Instructions

1. Clone the repository:
   ```bash
   git clone https://github.com/your-repository/synergy-simulation.git
   cd synergy-simulation
   ```

2. Open the desired RMarkdown file (`synergy_analysis_pipeline.Rmd` or others) in RStudio.

3. Execute the code:
   - By default, many chunks have `eval = FALSE`. To run specific sections, enable `eval = TRUE` in the desired code chunks.

4. Modify parameters as needed:
   - Adjust parameters like `nsyn`, `nrun`, or `param` to suit your simulation requirements.

5. Save outputs:
   - Results will be saved in the `RoMMS/simulation/sim_result/` directory.

---

## File Descriptions

### Main Files

- **`data_anslysis.html`**:
  Core pipeline for computing individual and population synergies and visualizing results with actual data.

- **`simulation_samplesize.html`**:
  Generates synthetic muscle activity data with varying sample sizes and computes synergies.

- **`simulation_outlier.html`**
  Generates synthetic muscle activity data with outliers and computes synergies.
  
### Scripts

- **`simulation_function.R`**:
  Contains utility functions for generating simulated muscle activity.

- **`modules.R`**:
  Additional helper functions for synergy computation and matching.

### Outputs

- Results are saved in `.rda` format under:
  - `results/`

- Visualizations are saved as `.pdf` files under:
  - `fig/`

---

## Visualization Outputs

1. **Synergy Comparisons**:
   - Visualizes `W` (activation patterns) and `H` (muscle contributions) matrices across methods (`Individual`, `Population`, `Concatenate`).

2. **Reconstruction Errors**:
   - Compares residuals of reconstructed synergies using bar plots.

3. **Estimation Errors**:
   - Evaluates errors in `W` and `H` matrices across sample sizes using smooth plots.
