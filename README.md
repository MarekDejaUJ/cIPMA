# cIPMA: Combined Importance-Performance Map Analysis

[![R-CMD-check](https://img.shields.io/badge/R--CMD--check-passing-brightgreen)](https://github.com/MarekDejaUJ/cIPMA)
[![License: GPL-3](https://img.shields.io/badge/License-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

The **cIPMA** package implements the Combined Importance-Performance Map Analysis framework proposed by [Hauff et al. (2024)](https://doi.org/10.1016/j.jretconser.2024.103723). It integrates PLS-SEM importance-performance analysis (via `seminr`) with Necessary Condition Analysis (via `NCA`) to jointly assess sufficiency and necessity of predictors in structural equation models.

## Features

- **Automated workflow**: Indicator rescaling, bottleneck computation, and visualization in a single function
- **Assumption diagnostics**: `check_assumptions()` validates reliability, HTMT, VIF, and sample size
- **Publication-ready plots**: cIPMA bubble plots with customizable aesthetics
- **Full output access**: PLS model, NCA results, bottleneck tables, and rescaled scores
- **Experimental structure learning**: Optional `discovery_cipma()` for exploratory analysis

## Installation

```r
# Install from GitHub
# install.packages("devtools")
devtools::install_github("MarekDejaUJ/cIPMA")
```

## Quick Start

```r
library(cIPMA)
library(seminr)

# 1. Load data
tam_data <- read.csv("tam.csv", check.names = FALSE)

# 2. Define measurement model
tam_mm <- constructs(
  composite("Perceived_Usefulness", multi_items("PU_0", 1:3)),
  composite("Compatibility",        multi_items("CO_0", 1:3)),
  composite("Ease_of_Use",          multi_items("EOU_0", 1:3)),
  composite("Emotional_Value",      multi_items("EMV_0", 1:3)),
  composite("Adoption_Intention",   multi_items("AD_0", 1:3)),
  composite("Technology_Use",       single_item("USE_01"))
)

# 3. Define structural model (theory-driven)
tam_sm <- relationships(
  paths(from = c("Compatibility", "Perceived_Usefulness", "Ease_of_Use", "Emotional_Value"), 
        to = "Adoption_Intention"),
  paths(from = c("Adoption_Intention", "Compatibility", "Perceived_Usefulness", "Ease_of_Use", "Emotional_Value"), 
        to = "Technology_Use")
)

# 4. Estimate and bootstrap PLS-SEM model
tam_pls <- estimate_pls(tam_data, tam_mm, tam_sm)
boot_tam <- bootstrap_model(tam_pls, nboot = 1000)

# 5. Define indicator scales (theoretical min/max)
my_scales <- list(
  PU_01 = c(1, 5), PU_02 = c(1, 5), PU_03 = c(1, 5),
  CO_01 = c(1, 5), CO_02 = c(1, 5), CO_03 = c(1, 5),
  EOU_01 = c(1, 5), EOU_02 = c(1, 5), EOU_03 = c(1, 5),
  EMV_01 = c(1, 5), EMV_02 = c(1, 5), EMV_03 = c(1, 5),
  AD_01 = c(1, 5), AD_02 = c(1, 5), AD_03 = c(1, 5),
  USE_01 = c(1, 7)
)

# 6. Run cIPMA
result <- cipma(
  model = boot_tam,
  target_construct = "Technology_Use",
  data = tam_data,
  scales = my_scales,
  target_level = 85,
  pls_model = tam_pls
)

# 7. View results
print(result)
plot(result)
check_assumptions(result)
```

## Output Components

```r
result$cIPMA_results      # Combined PLS-SEM and NCA statistics (Table A)
result$Bottleneck_summary # Bottleneck at target level (Table B)
result$Bottleneck_full    # Full bottleneck table (0-100% in 5% steps)
result$NCA_summary        # NCA effect sizes and p-values
result$Data_Rescaled      # Construct scores rescaled to 0-100
result$PLS_Model          # Base seminr model (for diagnostics)
result$Boot_Model         # Bootstrapped seminr model
```

## Assumption Checking

```r
check_assumptions(result)

# Output includes:
# 1. Positive outer weights
# 2. Reliability (alpha, rhoC, AVE, rhoA)
# 3. Discriminant validity (HTMT)
# 4. Multicollinearity (VIF)
# 5. Sample size adequacy
```

---

## ⚠️ Experimental: Structure Learning

> **Warning**: The `discovery_cipma()` function is **experimental** and intended for **exploratory analysis only**. Results should be interpreted with caution.

### Methodological Limitations

The `discovery_cipma()` function provides an interface to structure learning algorithms from `bnlearn`. However, users should be aware of important limitations:

1. **Exploratory only**: Discovered structures represent statistical associations, not proven causal relationships. They require theoretical justification before interpretation.

2. **Direct effects focus**: The current implementation analyzes only **direct predictors** of the target construct. Importance scores use PLS-SEM **total effects** (direct + indirect), while NCA operates **bivariately** between each predictor and outcome under ceteris paribus assumptions.

3. **No cascaded necessity**: NCA bottlenecks are not propagated through mediated causal chains. Whether sequential NCA along discovered pathways would yield different conclusions remains an open methodological question.

4. **Edge orientation**: Structure learning algorithms identify Markov Equivalence Classes, not unique DAGs. Users must inject background knowledge via `whitelist` and `blacklist` to resolve edge directions.

5. **Probabilistic vs. necessity logic**: Structure learning uses probabilistic scores (BIC) based on conditional independence, which may exclude variables that are statistically "unimportant" but critical as necessary conditions.

### Available Algorithms

`discovery_cipma()` supports 12 algorithms:

| Type | Algorithms |
|------|------------|
| **Score-based** | `hc` (Hill-Climbing), `tabu` (Tabu Search) |
| **Constraint-based** | `pc`, `gs`, `iamb`, `fast.iamb`, `inter.iamb`, `mmpc`, `hpc` |
| **Hybrid** | `mmhc` (recommended), `h2pc`, `rsmax2` |

### Example Usage

```r
library(bnlearn)  # Required for discovery_cipma

# Define constraints based on theory
blacklist <- data.frame(
  from = rep("Technology_Use", 5),
  to = c("Perceived_Usefulness", "Compatibility", "Ease_of_Use", 
         "Emotional_Value", "Adoption_Intention")
)

whitelist <- data.frame(
  from = c("Perceived_Usefulness", "Ease_of_Use", "Adoption_Intention"),
  to = c("Adoption_Intention", "Perceived_Usefulness", "Technology_Use")
)

# Run with hybrid algorithm
result_discovery <- discovery_cipma(
  data = tam_data,
  measurement_model = tam_mm,
  target_construct = "Technology_Use",
  scales = my_scales,
  target_level = 85,
  algorithm = "mmhc",
  blacklist = blacklist,
  whitelist = whitelist,
  nboot = 500
)

# Run with constraint-based algorithm (Optimized)
result_gs <- discovery_cipma(
    data = tam_data,
    measurement_model = tam_mm,
    target_construct = "Technology_Use",
    scales = my_scales,
    target_level = 80,
    algorithm = "gs",         # The slow algorithm
    blacklist = blacklist,
    whitelist = whitelist,
    nboot = 500,
    algorithm.args = list(
        max.sx = 5,
        alpha = 0.01
    )
)

# View discovered structure
print(result_discovery$discovered_arcs)
plot(result_discovery$learned_graph)
print(result_discovery)

print(result_gs$discovered_arcs)
plot(result_gs$learned_graph)
print(result_gs)
```
### Performance Note (Constraint-Based Algorithms)
Algorithms like gs (Grow-Shrink) or iamb rely on conditional independence tests. When applied to dense PLS-SEM datasets where constructs are highly correlated, these algorithms may become computationally prohibitive (taking hours or causing crashes).
Recommended optimization for 'gs', 'iamb', or 'pc': Use algorithm.args to limit the size of the conditioning set (max.sx). This forces the algorithm to approximate, significantly improving speed

### Recommendation

For rigorous analysis, use the theory-driven `cipma()` function with a pre-specified structural model. Reserve `discovery_cipma()` for:
- Exploratory research where theory is underdeveloped
- Comparing theory-driven vs. data-driven structures
- Identifying potentially missing paths for further theoretical consideration

---

## References

- Hauff, S., Richter, N. F., Sarstedt, M., & Ringle, C. M. (2024). Importance and performance in PLS-SEM and NCA: Introducing the combined importance-performance map analysis (cIPMA). *Journal of Retailing and Consumer Services*, 78, 103723.

- Dul, J. (2016). Necessary Condition Analysis (NCA): Logic and methodology of "necessary but not sufficient" causality. *Organizational Research Methods*, 19(1), 10-52.

- Dul, J. (2021). *Conducting Necessary Condition Analysis*. Sage Publications.

- Ray, S., Danks, N., & Calero Valdez, A. (2021). SEMinR: Domain-specific language for building, estimating, and visualizing structural equation models in R. *SSRN Electronic Journal*.

- Scutari, M. (2010). Learning Bayesian networks with the bnlearn R package. *Journal of Statistical Software*, 35(3), 1-22.

## Citation

If you use this package, please cite:

```
Deja, M. (2025). cIPMA: Combined Importance-Performance Map Analysis. 
R package version 0.1.6. https://github.com/MarekDejaUJ/cIPMA
```

