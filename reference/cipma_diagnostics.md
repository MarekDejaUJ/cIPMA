# QNCA Diagnostics of a cIPMA

Runs the diagnostics of Quantile Necessary Condition Analysis on every
predictor of a cIPMA, with the rescaled scores and the data-range scope
of each pair. The reference band compares the observed effect size with
a Gaussian-copula benchmark that keeps the marginals and the rank
dependence of the pair and imposes no attainability floor. The deletion
screen reports the change of the effect size when one case, or `k` cases
jointly, leave the sample. The resolution table gives, target by target,
the conditioning-set size, the selected rank, the inherited flag and the
resistance of the frontier at the chosen tolerance. The consistency
probe follows the effect size across subsample sizes. The addition
screen, run when `points` is given, adds cases at declared positions.

## Usage

``` r
cipma_diagnostics(
  x,
  constructs = NULL,
  pi_grid = c(1, 0.95, 0.9),
  k = 1L,
  M = 199L,
  probe = TRUE,
  reps = 50L,
  points = NULL,
  counts = 0:5,
  n_grid = 50L,
  seed = NULL
)
```

## Arguments

- x:

  A cipma object.

- constructs:

  Predictors to diagnose; default all.

- pi_grid:

  Tolerances of the band, the deletion screen and the addition screen.
  Default 1, 0.95, 0.90.

- k:

  Number of cases deleted jointly in the deletion screen. Default 1.

- M:

  Number of benchmark datasets of the reference band. Default 199.

- probe:

  Logical. Run the consistency probe. Default `TRUE`.

- reps:

  Subsamples per size of the consistency probe. Default 50.

- points:

  Optional `c(x, y)` position, or two-column matrix of positions, of
  cases added by the addition screen.

- counts:

  Numbers of added cases of the addition screen. Default 0 to 5.

- n_grid:

  Number of outcome levels. Default 50.

- seed:

  Optional seed; the random number state of the session is restored
  afterwards.

## Value

An object of class `cipma_diagnostics`: per construct, the reference
band, the deletion screen, the resolution table, the consistency probe
and the addition screen, with the settings.

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
res <- cipma(boot, "Satisfaction", mobi, scales, target_level = 80, B = 99, seed = 1)
diag <- cipma_diagnostics(res, M = 49, reps = 10, seed = 1)
diag
#> QNCA diagnostics of the cIPMA predictors of Satisfaction (tolerances 1.00, 0.95, 0.90)
#> 
#> Reference band (49 benchmark datasets): observed d, upper band limit, excess
#>  Construct   pi d_obs upper excess p_null Above_band
#>      Image 1.00 0.297 0.378 -0.080  0.660         no
#>      Image 0.95 0.476 0.526 -0.050  0.700         no
#>      Image 0.90 0.553 0.568 -0.015  0.380         no
#>      Value 1.00 0.000 0.251 -0.251  1.000         no
#>      Value 0.95 0.381 0.460 -0.079  0.940         no
#>      Value 0.90 0.462 0.519 -0.057  1.000         no
#> 
#> Deletion screen (k = 1): most influential deletion and its relative change of d per tolerance
#>  Construct Cases       Relative_change Flagged_rows
#>      Image   242 0.067 / 0.016 / 0.003           13
#>      Value    73    NA / 0.016 / 0.018           10
#> 
#> Resolution at pi = 1.00 for outcome targets up to 80
#>  Construct Targets Inherited Smallest_resistance Smallest_attaining_set
#>      Image      36         0                   0                     79
#>      Value      36         0                   0                     79
#> 
#> Consistency probe: divergence of the drift slopes and the slopes at the two tolerances
#>  Construct Divergence              Drift
#>      Image      7.582  7.87234 / 0.28995
#>      Value     18.702 20.62543 / 1.92349
#> 
#> An empty zone above the band exceeds the dependence-matched benchmark; the band is a
#> proposed diagnostic with operating characteristics still to be established.
# }
```
