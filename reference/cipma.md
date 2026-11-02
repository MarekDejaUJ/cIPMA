# Run Combined Importance-Performance Map Analysis (cIPMA)

Combines the importance-performance map analysis of a PLS-SEM model with
Quantile Necessary Condition Analysis (QNCA) of each predictor of the
target construct (Hauff et al. 2024). Importance is the total effect on
the target construct, and performance is the mean of the latent variable
scores rescaled to 0-100. Necessity is decided at the strict limit
`pi = 1`, where the QNCA frontier equals the CE-FDH ceiling of Necessary
Condition Analysis: a predictor is necessary when its CE-FDH effect size
reaches 0.10 and its permutation p-value stays below 0.05. The tolerance
`pi` sets the share `1 - pi` of cases that attain the target level and
may sit below the reported requirement; the bottleneck, the share of
cases failing it and the effect size at `pi` follow the chosen
tolerance. The tolerance path repeats the effect size and the bottleneck
for every tolerance in `pi_path`.

## Usage

``` r
cipma(
  model,
  target_construct,
  data,
  scales,
  target_level = 85,
  nca_rep = NULL,
  pls_model = NULL,
  pi = 1,
  pi_path = seq(1, 0.9, by = -0.01),
  B = 1999L,
  n_grid = 400L,
  n_grid_test = 100L,
  seed = NULL
)
```

## Arguments

- model:

  A bootstrapped seminr model object.

- target_construct:

  String. The name of the target construct.

- data:

  The raw data frame used for the model.

- scales:

  A named list of theoretical minimum and maximum per indicator.

- target_level:

  Numeric. The desired level of the target construct on the 0-100 scale.
  Default 85.

- nca_rep:

  Deprecated. Number of permutations; use `B`.

- pls_model:

  Optional. The base (non-bootstrapped) seminr model for the assumption
  checks.

- pi:

  Tolerance in (0, 1\] for the bottleneck and the effect size at the
  chosen tolerance. Default 1, which reproduces the CE-FDH results of
  NCA.

- pi_path:

  Tolerances of the tolerance path. Default 1.00 to 0.90 in steps of
  0.01.

- B:

  Number of permutations of the significance test. Default 1999; 0 skips
  the test.

- n_grid:

  Number of outcome levels of the QNCA effect sizes. Default 400.

- n_grid_test:

  Number of outcome levels inside the permutation test. Default 100.

- seed:

  Optional seed of the permutation test; the random number state of the
  session is restored afterwards.

## Value

An object of class `cipma`: a list with the combined results
(`cIPMA_results`), the bottleneck at the target level
(`Bottleneck_summary`), the bottleneck table across outcome levels
(`Bottleneck_full`), the tolerance path (`Tolerance_Path`), the
settings, the rescaled scores and the models.

## References

Hauff, S., Richter, N. F., Sarstedt, M., & Ringle, C. M. (2024).
Importance and performance in PLS-SEM and NCA: Introducing the combined
importance-performance map analysis (cIPMA). Journal of Retailing and
Consumer Services, 78, 103723.
[doi:10.1016/j.jretconser.2024.103723](https://doi.org/10.1016/j.jretconser.2024.103723)

## Examples

``` r
# \donttest{
library(seminr)
mm <- constructs(
  composite("Image", multi_items("IMAG", 1:5)),
  composite("Expectation", multi_items("CUEX", 1:3)),
  composite("Value", multi_items("PERV", 1:2)),
  composite("Satisfaction", multi_items("CUSA", 1:3))
)
sm <- relationships(
  paths(from = c("Image", "Expectation", "Value"), to = "Satisfaction")
)
pls <- estimate_pls(mobi, mm, sm)
#> Generating the seminr model
#> All 250 observations are valid.
boot <- bootstrap_model(pls, nboot = 50, cores = 1, seed = 1)
#> Bootstrapping model using seminr...
#> SEMinR Model successfully bootstrapped
scales <- setNames(rep(list(c(1, 10)), ncol(mobi)), names(mobi))
res <- cipma(boot, "Satisfaction", mobi, scales, target_level = 80,
             pi = 0.95, B = 99, seed = 1, pls_model = pls)
res
#> cIPMA: target Satisfaction at level 80 (n = 250, tolerance pi = 0.95)
#> 
#>    Construct Importance Performance     d      p Necessary Required Failing
#>        Image      0.447        73.6 0.295 0.0100       yes   64.334   24.0%
#>  Expectation      0.172        72.5 0.193 0.0100       yes   58.892   14.4%
#>        Value      0.316        62.4 0.000 1.0000        no   44.444   14.4%
#> 
#> Necessary: CE-FDH d >= 0.10 and p < 0.05 at pi = 1.00. Required and Failing:
#> bottleneck at pi = 0.95 (NN: no requirement; NA: level not attained).
#> summary(x) adds the readings; summary(x, detail = "full") adds the tolerance path and
#> the bottleneck table; cipma_diagnostics(x) runs the QNCA diagnostics.
summary(res)
#> Combined importance-performance map analysis
#> Target construct: Satisfaction, target level: 80, cases: 250
#> Tolerance of the bottleneck: pi = 0.95; permutations: 99; outcome levels: 400 (test: 100); qnca 0.4.0
#> 
#> Combined map
#>    Construct Importance Performance     d      p Necessary Required Failing
#>        Image      0.447        73.6 0.295 0.0100       yes   64.334   24.0%
#>  Expectation      0.172        72.5 0.193 0.0100       yes   58.892   14.4%
#>        Value      0.316        62.4 0.000 1.0000        no   44.444   14.4%
#> 
#> PLS-SEM significance, effect at the chosen tolerance, failing cases
#>    Construct   PLS_p d_0.95 p_0.95 Failing_cases
#>        Image <0.0001  0.475 0.0100            60
#>  Expectation  0.0076  0.400 0.0100            36
#>        Value <0.0001  0.381 0.0100            36
#> 
#> Readings
#> - Image is necessary at the strict limit (d = 0.295, p = 0.0100); the required level of Image for Satisfaction >= 80 at pi = 0.95 is 64.334, with 24.0% of the cases below the required level.
#> - Expectation is necessary at the strict limit (d = 0.193, p = 0.0100); the required level of Expectation for Satisfaction >= 80 at pi = 0.95 is 58.892, with 14.4% of the cases below the required level.
#> - Value is not necessary at the strict limit (d = 0.000, p = 1.0000); the required level of Value for Satisfaction >= 80 at pi = 0.95 is 44.444, with 14.4% of the cases below the required level.
#> 
#> Effect sizes at pi = 0.95 are relative to the tolerance and rise as pi falls;
#> the NCA landmark of 0.10 applies to the CE-FDH effect size at pi = 1.00.
# }
```
