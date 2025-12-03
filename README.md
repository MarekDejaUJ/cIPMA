# cIPMA: Combined Importance-Performance Map Analysis

The **cIPMA** package implements the Combined Importance-Performance Map Analysis approach proposed by [Hauff et al. (2024)](https://doi.org/10.1016/j.jretconser.2024.103723). It integrates results from PLS-SEM (using `seminr`) with Necessary Condition Analysis (using `NCA`) to identify both sufficiency and necessity factors.

## Installation

You can install the development version of cIPMA like so:

```r
# install.packages("devtools")
devtools::install_github("yourusername/cIPMA")
```

## Example Usage

Here is a full workflow using the extended TAM dataset (included as `tam.csv` in your project).

```r
library(cIPMA)
library(seminr)
library(NCA)
library(ggplot2)
library(ggrepel)

# 1. Load Data
# Ensure "tam.csv" is in your working directory
tam_data <- read.csv("tam.csv", check.names = FALSE)

# 2. Define Seminr Model
tam_mm <- constructs(
  composite("Perceived_Usefulness", multi_items("PU_0", 1:3)),
  composite("Compatibility",        multi_items("CO_0", 1:3)),
  composite("Ease_of_Use",          multi_items("EOU_0", 1:3)),
  composite("Emotional_Value",      multi_items("EMV_0", 1:3)),
  composite("Adoption_Intention",   multi_items("AD_0", 1:3)),
  composite("Technology_Use",       single_item("USE_01"))
)

tam_sm <- relationships(
  paths(from = c("Compatibility", "Perceived_Usefulness", "Ease_of_Use", "Emotional_Value"), to = "Adoption_Intention"),
  paths(from = c("Adoption_Intention", "Compatibility", "Perceived_Usefulness", "Ease_of_Use", "Emotional_Value"), to = "Technology_Use")
)

# 3. Estimate & Bootstrap (REQUIRED for cIPMA CIs)
# Note: cIPMA requires a bootstrapped model for p-values
tam_pls <- estimate_pls(tam_data, tam_mm, tam_sm)
boot_tam <- bootstrap_model(tam_pls, nboot = 1000)

# 4. Define Scales for Rescaling (Theoretical Min/Max)
my_scales <- list(
  PU_01 = c(1, 5), PU_02 = c(1, 5), PU_03 = c(1, 5),
  CO_01 = c(1, 5), CO_02 = c(1, 5), CO_03 = c(1, 5),
  EOU_01 = c(1, 5), EOU_02 = c(1, 5), EOU_03 = c(1, 5),
  EMV_01 = c(1, 5), EMV_02 = c(1, 5), EMV_03 = c(1, 5),
  AD_01 = c(1, 5), AD_02 = c(1, 5), AD_03 = c(1, 5),
  USE_01 = c(1, 7)
)

# 5. Run cIPMA
result <- cipma(
  model = boot_tam, 
  target_construct = "Technology_Use", 
  data = tam_data, 
  scales = my_scales, 
  target_level = 85
)

# 6. Outputs
print(result) # Prints Table A (Stats) and Table B (Bottlenecks)
my_cipma_plot <- plot(result) # Generates the Figure

# 7. Save Plot
ggsave(
  filename = "cipma_plot.png",
  plot = my_cipma_plot,
  width = 700, 
  height = 550, 
  units = "px",
  dpi = 100,
  bg = "white"
)
```

### Data-Driven Analysis (Causal Discovery)

If you do not have a pre-defined structural model, you can use the `discovery_cipma` function. [cite_start]This uses the Hill-Climbing algorithm from the `bnlearn` package [cite: 11, 3316] to learn the causal structure from your data before running the standard cIPMA.

```r
# Run Data-Driven cIPMA
result2 <- discovery_cipma(
  data              = tam_data,
  measurement_model = tam_mm,
  target_construct  = "Technology_Use",
  scales            = my_scales,  # your PU/CO/EOU/EMV/AD/USE scales
  target_level      = 85
)
print(result2)

plot(result2)
```

## References
* Hauff, S., Richter, N. F., Sarstedt, M., & Ringle, C. M. (2024). Importance and performance in PLS-SEM and NCA: Introducing the combined importance-performance map analysis (cIPMA). *Journal of Retailing and Consumer Services*, 78, 103723.
* Dul, J. (2016). Necessary Condition Analysis (NCA). Logic and Methodology of 'Necessary but not Sufficient' causality. *Organizational Research Methods*, 19(1), 10-52.
