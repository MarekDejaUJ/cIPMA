# Summarise cIPMA Results

Collects the results of a cIPMA at three levels of detail. `"brief"`
gives the combined table. `"standard"` adds the PLS-SEM p-values, the
effect size and p-value at the chosen tolerance and one reading per
construct. `"full"` adds the tolerance path of the effect size and of
the required level, the bottleneck table across outcome levels and the
settings.

## Usage

``` r
# S3 method for class 'cipma'
summary(object, detail = c("standard", "brief", "full"), digits = 3, ...)
```

## Arguments

- object:

  A cipma object.

- detail:

  One of `"standard"`, `"brief"` or `"full"`.

- digits:

  Decimals of effect sizes and importance.

- ...:

  Ignored.

## Value

An object of class `summary.cipma`, printed by its print method.

## Examples

``` r
# \donttest{
library(seminr)
mm <- constructs(
  composite("Image", multi_items("IMAG", 1:5)),
  composite("Value", multi_items("PERV", 1:2)),
  composite("Satisfaction", multi_items("CUSA", 1:3))
)
sm <- relationships(paths(from = c("Image", "Value"), to = "Satisfaction"))
pls <- estimate_pls(mobi, mm, sm)
#> Generating the seminr model
#> All 250 observations are valid.
boot <- bootstrap_model(pls, nboot = 50, cores = 1, seed = 1)
#> Bootstrapping model using seminr...
#> SEMinR Model successfully bootstrapped
scales <- setNames(rep(list(c(1, 10)), ncol(mobi)), names(mobi))
res <- cipma(boot, "Satisfaction", mobi, scales, target_level = 80,
             pi = 0.95, B = 99, seed = 1)
summary(res, detail = "full")
#> Combined importance-performance map analysis
#> Target construct: Satisfaction, target level: 80, cases: 250
#> Tolerance of the bottleneck: pi = 0.95; permutations: 99; outcome levels: 400 (test: 100); qnca 0.4.0
#> 
#> Combined map
#>  Construct Importance Performance     d      p Necessary Required Failing
#>      Image      0.521        73.6 0.296 0.0100       yes   64.343   24.0%
#>      Value      0.342        62.4 0.000 1.0000        no   44.444   14.4%
#> 
#> PLS-SEM significance, effect at the chosen tolerance, failing cases
#>  Construct   PLS_p d_0.95 p_0.95 Failing_cases
#>      Image <0.0001  0.475 0.0100            60
#>      Value <0.0001  0.381 0.0100            36
#> 
#> Readings
#> - Image is necessary at the strict limit (d = 0.296, p = 0.0100); the required level of Image for Satisfaction >= 80 at pi = 0.95 is 64.343, with 24.0% of the cases below the required level.
#> - Value is not necessary at the strict limit (d = 0.000, p = 1.0000); the required level of Value for Satisfaction >= 80 at pi = 0.95 is 44.444, with 14.4% of the cases below the required level.
#> 
#> Effect sizes at pi = 0.95 are relative to the tolerance and rise as pi falls;
#> the NCA landmark of 0.10 applies to the CE-FDH effect size at pi = 1.00.
#> 
#> Tolerance path: effect size d
#>    pi Image Value
#>  1.00 0.296 0.000
#>  0.99 0.351 0.178
#>  0.98 0.401 0.262
#>  0.97 0.438 0.337
#>  0.96 0.460 0.364
#>  0.95 0.475 0.381
#>  0.94 0.491 0.393
#>  0.93 0.517 0.403
#>  0.92 0.527 0.437
#>  0.91 0.545 0.450
#>  0.90 0.553 0.462
#> 
#> Tolerance path: required level for Satisfaction >= 80
#>    pi  Image  Value
#>  1.00 57.807     NN
#>  0.99 57.807 26.714
#>  0.98 60.851 33.333
#>  0.97 61.924 40.012
#>  0.96 63.065 44.444
#>  0.95 64.343 44.444
#>  0.94 65.825 44.444
#>  0.93 65.825 44.503
#>  0.92 66.667 51.123
#>  0.91 69.103 51.182
#>  0.90 72.579 53.369
#> 
#> Bottleneck table at pi = 0.95 (NN: no requirement; NA: level not attained)
#>  Outcome_Level  Image  Value
#>              0 49.188 26.714
#>             10 49.188 26.714
#>             20 49.188 26.714
#>             30 49.329 26.714
#>             40 49.896 31.146
#>             50 54.167 31.146
#>             60 59.026 37.825
#>             70 61.677 44.444
#>             80 64.343 44.444
#>             90 79.370 44.444
#>            100 81.420 44.444
# }
```
